#!/usr/bin/env Rscript
# Fixed-parameter replication of GSM6258256 on GSM6258257 and GSM6258258.
options(width=180, future.globals.maxSize=8*1024^3)
.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5','G:/gurobi/Rlib45',.libPaths()))
args <- commandArgs(trailingOnly=TRUE)
sample_id <- if(length(args)) args[1] else stop('Provide GSM accession')
stopifnot(sample_id %in% c('GSM6258257','GSM6258258'))
root <- 'G:/1a空转rperio'
out <- file.path(root,sample_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
setwd(out)
pkgs <- c('Seurat','Matrix','spacexr','SpatialExperiment','SummarizedExperiment','mistyR','ggplot2','patchwork','R.utils')
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)))
write.csv(data.frame(package=pkgs,version=vapply(pkgs,function(p) as.character(packageVersion(p)),character(1))), 'package_versions.csv',row.names=FALSE)
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(patchwork)})
set.seed(20260824)
save_plot <- function(name,p,w=7,h=6) {
  ggsave(paste0(name,'.pdf'),p,width=w,height=h)
  ggsave(paste0(name,'.png'),p,width=w,height=h,dpi=220)
}
rawfiles <- list.files(file.path(root,'GSE206621_RAW'),pattern=paste0('^',sample_id,'_'),full.names=TRUE)
stopifnot(length(rawfiles)>=7)
write.csv(data.frame(file=basename(rawfiles),bytes=file.info(rawfiles)$size,md5=unname(tools::md5sum(rawfiles))), 'input_manifest.csv',row.names=FALSE)
dir.create('counts',showWarnings=FALSE);dir.create('spatial',showWarnings=FALSE)
for(f in rawfiles) {
  nm <- sub(paste0('^',sample_id,'_[^_]+_'),'',basename(f))
  target <- file.path(if(nm %in% c('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz')) 'counts' else 'spatial',nm)
  if(grepl('^spatial',target)) {
    dest <- sub('[.]gz$','',target)
    if(!file.exists(dest)) R.utils::gunzip(f,destname=dest,remove=FALSE,overwrite=FALSE)
  } else if(!file.exists(target)) stopifnot(file.copy(f,target,overwrite=FALSE))
}
if(!file.exists('spatial_PD_input.rda')) {
  counts <- Read10X('counts')
  if(is.list(counts)) counts <- counts[['Gene Expression']]
  positions <- read.csv('spatial/tissue_positions_list.csv',header=FALSE,stringsAsFactors=FALSE)
  colnames(positions) <- c('spot_id','in_tissue','array_row','array_col','pxl_row','pxl_col')
  tissue <- intersect(colnames(counts),positions$spot_id[positions$in_tissue==1])
  stopifnot(length(tissue)>100,!anyDuplicated(positions$spot_id))
  n_matrix <- ncol(counts)
  stRNA <- CreateSeuratObject(counts=counts[,tissue,drop=FALSE],assay='Spatial',project=sample_id)
  rm(counts);gc()
  stRNA[['slice1']] <- Read10X_Image('spatial',filter.matrix=TRUE)
  stRNA$sample_id <- sample_id;stRNA$condition <- 'PD'
  stRNA[['percent.mt']] <- PercentageFeatureSet(stRNA,pattern='^MT-',assay='Spatial')
  write.csv(stRNA[[]],'spot_QC_metrics.csv')
  write.csv(data.frame(sample=sample_id,stage=c('matrix_barcodes','in_tissue'),n_spots=c(n_matrix,ncol(stRNA))), 'QC_spot_counts.csv',row.names=FALSE)
  save_plot('QC_metrics',VlnPlot(stRNA,features=c('nCount_Spatial','nFeature_Spatial','percent.mt'),pt.size=0.1,ncol=3),12,4.5)
  # Match the original PD workflow: in-tissue filter; no new count/mito thresholds.
  stRNA <- SCTransform(stRNA,assay='Spatial',verbose=FALSE,seed.use=20260824)
  stRNA <- RunPCA(stRNA,assay='SCT',verbose=FALSE,seed.use=20260824)
  stRNA <- FindNeighbors(stRNA,dims=1:15,verbose=FALSE)
  stRNA <- FindClusters(stRNA,resolution=0.8,algorithm=1,random.seed=20260824,verbose=FALSE)
  stRNA <- RunUMAP(stRNA,dims=1:15,seed.use=20260824,verbose=FALSE)
  # These are unsupervised spatial clusters, NOT transferred anatomical regions.
  stRNA$Region <- paste0('Spatial_cluster_',stRNA$seurat_clusters)
  save(stRNA,file='spatial_PD_input.rda')
  save_plot('spatial_clusters',SpatialDimPlot(stRNA,label=TRUE,label.size=4,pt.size.factor=2)+ggtitle(sample_id),7,6)
  save_plot('UMAP_clusters',DimPlot(stRNA,label=TRUE,pt.size=1.5)+theme_classic(base_size=14)+theme(aspect.ratio=1,panel.border=element_rect(fill=NA,color='black')),7,6)
  DefaultAssay(stRNA) <- 'Spatial'
  stRNA <- NormalizeData(stRNA,verbose=FALSE)
  genes <- intersect(c('PI16','AGT','C3','POSTN','COL1A1','FN1','TGFB1','CD44'),rownames(stRNA))
  save_plot('mechanism_gene_spatial',SpatialFeaturePlot(stRNA,features=genes,ncol=4,pt.size.factor=2),16,8)
  expr <- FetchData(stRNA,vars=genes,layer='data')
  write.csv(data.frame(spot_id=rownames(expr),expr,check.names=FALSE),'mechanism_gene_log_normalized.csv',row.names=FALSE)
  write.csv(data.frame(gene=genes,n_spots=nrow(expr),positive_percent=colMeans(expr>0)*100,mean_log_normalized=colMeans(expr)),'mechanism_gene_summary.csv',row.names=FALSE)
  saveRDS(stRNA,'spatial_normalized.rds')
  rm(stRNA);gc()
}
# Reuse exact original reference selection, random cap, createRctd/runRctd settings.
# Stop before region tests: new cluster IDs are not anatomical annotations.
if(!file.exists('RCTD_PD_refinedFib_ECM_results.rds')) {
  src <- readLines(file.path(out,'provenance','RCTD_PD_refinedFib_ECM_rerun.R'),warn=FALSE,encoding='UTF-8')
  cutoff <- which(grepl('^region_plot <-',src))[1]-1L
  run_rctd <- function() eval(parse(text=src[seq_len(cutoff)]),envir=environment())
  run_rctd()
  while(sink.number()>0) sink()
  gc()
}
if(!file.exists('MISTy_PD_RCTD_refinedFib_ECM_collect.rds')) {
  run_misty_original <- function() sys.source(file.path(out,'provenance','MISTy_PD_refinedFib_ECM_full.R'),envir=environment())
  run_misty_original()
  while(sink.number()>0) sink()
}
# Explicit coordinate/neighborhood audit for the legacy fixed-distance MISTy model.
st <- readRDS('stRNA_RCTD_refinedFib_ECM_PD.rds')
ca <- attributes(attributes(attr(st,'images')[[1]])$boundaries[['centroids']])
xy <- as.matrix(ca$coords);rownames(xy) <- ca$cells
w <- read.csv('RCTD_PD_refinedFib_ECM_weights_full.csv',check.names=FALSE)
w <- w[rowSums(w[,-1,drop=FALSE])>0,,drop=FALSE]
xy <- xy[w$spot_id,,drop=FALSE]
dd <- as.matrix(dist(xy));diag(dd) <- Inf
write.csv(data.frame(sample=sample_id,n_spots=nrow(xy),median_nearest_distance=median(apply(dd,1,min)),mean_neighbors_50=mean(rowSums(dd<50)),isolated_fraction_50=mean(rowSums(dd<50)==0),note='Distances use original image-centroid coordinates; spots are not biological replicates'),'coordinate_neighborhood_audit.csv',row.names=FALSE)
write.csv(data.frame(sample=sample_id,n_RCTD_output=nrow(read.csv('RCTD_PD_refinedFib_ECM_weights_full.csv')),n_positive_full_weights=nrow(w)),'RCTD_spot_counts.csv',row.names=FALSE)
capture.output(sessionInfo(),file='sessionInfo.txt')
writeLines(c('Analysis complete.','Primary replication uses the original GSM6258256 settings and the same refined single-cell reference.','No disease-group variable used for clustering. No new anatomical region labels assigned.','RNA maps: log-normalized raw Spatial counts. RCTD: integer raw counts.','RCTD focal maps use doublet weights; MISTy/correlation use continuous weights_full, as in the original analysis.','Correlations and their P values are descriptive spot-level calculations, NOT donor-level evidence; spatial autocorrelation is not adjusted.','Two additional PD sections provide within-dataset replication; donor independence remains to be verified.','Negative or discordant results must be retained. Co-localization/MISTy cannot establish direct signaling or causality.'),'ANALYSIS_NOTES.txt')
cat('COMPLETE:',sample_id,'\n')
