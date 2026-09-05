options(stringsAsFactors=FALSE)
user_library <- 'C:/Users/32266/AppData/Local/R/win-library/4.5'
.libPaths(c(user_library,'G:/gurobi/Rlib45',.libPaths()))
pkgs <- c('SeuratObject','Matrix','edgeR','ggplot2')
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Existing packages unavailable: ',paste(missing,collapse=', '),'. No installation attempted.')
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(edgeR);library(ggplot2)})
set.seed(20260902)
root <- 'G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final'
input <- file.path(root,'results/objects/fibroblast_revision_final_no_unassigned.rds')
out <- 'all_subset_model_audit_20260902/Figure6C_three_gene_GSM'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
figure_dir <- 'G:/1Yunvjian/1A\u7389\u5973\u714e/Figure'
targets <- c('AGT','CXCR4','FOS')
minimum_cells <- 50L

message('[1/4] Reading the revised Figure 6 fibroblast object')
fib <- readRDS(input)
stopifnot(inherits(fib,'Seurat'),ncol(fib)==16032L)
if(inherits(fib[['RNA']],'Assay5')) {
  count_layers <- Layers(fib[['RNA']],search='^counts')
  if(length(count_layers)>1L)fib[['RNA']]<-JoinLayers(fib[['RNA']],layers=count_layers,new='counts')
}
counts <- GetAssayData(fib,assay='RNA',layer='counts')
meta <- fib[[]]
stopifnot(identical(colnames(counts),rownames(meta)),all(c('sample','group','dataset','state_final')%in%names(meta)),all(targets%in%rownames(counts)))
stopifnot(!anyNA(meta[,c('sample','group','dataset','state_final')]),all(meta$group%in%c('NC','PD')))
stopifnot(setequal(unique(as.character(meta$state_final)),c('Inflammatory_Fib','Activated_Fib','ECM_Fib','Adventitial_Fib','PI16_Fib')))
if(inherits(counts,'sparseMatrix'))stopifnot(all(counts@x>=0),all(abs(counts@x-round(counts@x))<1e-8))
samples <- sort(unique(as.character(meta$sample)))
sm <- do.call(rbind,lapply(samples,function(s){
  z<-meta[as.character(meta$sample)==s,,drop=FALSE]
  stopifnot(length(unique(z$group))==1L,length(unique(z$dataset))==1L)
  data.frame(GSM=s,group=as.character(z$group[1]),dataset=as.character(z$dataset[1]),n_fibroblasts=nrow(z))
}))
sm$included <- sm$n_fibroblasts>=minimum_cells
sm$exclusion_reason <- ifelse(sm$included,'','Fewer than 50 qualified fibroblasts')
write.csv(sm,file.path(out,'Table_GSM_inclusion.csv'),row.names=FALSE)
write.csv(as.data.frame(table(sm$dataset,sm$group,sm$included)),file.path(out,'Table_dataset_group_inclusion.csv'),row.names=FALSE)

message('[2/4] Summing raw RNA counts within GSM and TMM normalizing')
membership <- sparseMatrix(i=seq_len(nrow(meta)),j=match(as.character(meta$sample),samples),x=1,dims=c(nrow(meta),length(samples)))
pb_all <- as.matrix(counts%*%membership)
rownames(pb_all)<-rownames(counts);colnames(pb_all)<-samples
stopifnot(sum(pb_all)==sum(counts),identical(sm$GSM,colnames(pb_all)))
eligible <- sm[sm$included,,drop=FALSE]
eligible$group <- factor(eligible$group,levels=c('NC','PD'))
eligible$dataset <- factor(eligible$dataset)
stopifnot(nrow(eligible)==28L,sum(eligible$group=='NC')==16L,sum(eligible$group=='PD')==12L)
pb <- pb_all[,eligible$GSM,drop=FALSE]
pb <- pb[rowSums(pb)>0,,drop=FALSE]
stopifnot(all(targets%in%rownames(pb)),all(colSums(pb)>0))
y <- DGEList(counts=pb)
y <- calcNormFactors(y,method='TMM')
logcpm <- cpm(y,log=TRUE,prior.count=2)
eligible$library_size <- y$samples$lib.size
eligible$TMM_factor <- y$samples$norm.factors
eligible$effective_library_size <- eligible$library_size*eligible$TMM_factor
plotdata <- do.call(rbind,lapply(targets,function(g){
  z<-eligible;z$gene<-g;z$raw_count<-as.numeric(pb[g,]);z$log2_CPM<-as.numeric(logcpm[g,]);z
}))
plotdata$gene <- factor(plotdata$gene,levels=targets)
stopifnot(nrow(plotdata)==84L,!anyDuplicated(paste(plotdata$GSM,plotdata$gene)),all(is.finite(plotdata$log2_CPM)))
write.csv(plotdata,file.path(out,'Table_Figure6C_GSM_expression.csv'),row.names=FALSE)
write.csv(eligible,file.path(out,'Table_TMM_normalization.csv'),row.names=FALSE)
summ <- do.call(rbind,lapply(split(plotdata,list(plotdata$gene,plotdata$group),drop=TRUE),function(z){
  data.frame(gene=as.character(z$gene[1]),group=as.character(z$group[1]),n_GSM=nrow(z),median=median(z$log2_CPM),Q1=unname(quantile(z$log2_CPM,.25)),Q3=unname(quantile(z$log2_CPM,.75)))
}))
write.csv(summ,file.path(out,'Table_descriptive_summary.csv'),row.names=FALSE)
saveRDS(list(counts=pb,sample_info=eligible,targets=targets,minimum_cells=minimum_cells),file.path(out,'GSM_pseudobulk_counts_and_metadata.rds'))
saveRDS(plotdata,file.path(out,'Figure6C_plot_data.rds'))

message('[3/4] Exporting a compact three-panel box-and-dot plot')
group_cols<-c(NC='#00BFC4',PD='#F8766D')
dataset_cols<-c(GSE152042='#E58B27',GSE164241='#00A582',GSE171213='#3B7DDD')
# Reuse exactly the same horizontal offset for each GSM in all three facets.
offset<-setNames(runif(nrow(eligible),-.12,.12),eligible$GSM)
plotdata$xpos<-as.integer(plotdata$group)+unname(offset[plotdata$GSM])
p<-ggplot(plotdata,aes(x=as.integer(group),y=log2_CPM))+
  geom_boxplot(aes(group=group,fill=group,colour=group),width=.54,linewidth=.65,outlier.shape=NA,alpha=.22,show.legend=FALSE)+
  geom_point(aes(x=xpos,colour=dataset),size=2.0,alpha=.9)+
  scale_fill_manual(values=group_cols)+
  scale_colour_manual(values=c(group_cols,dataset_cols),breaks=levels(eligible$dataset),name='Dataset')+
  scale_x_continuous(breaks=c(1,2),labels=c('NC','PD'),limits=c(.55,2.45))+
  scale_y_continuous(expand=expansion(mult=c(.07,.12)))+
  facet_wrap(~gene,nrow=1,scales='free_y')+
  labs(title='Three-gene expression in fibroblasts',subtitle='GSM-level pseudobulk | NC n = 16; PD n = 12',x=NULL,y='Expression (log2 CPM)')+
  theme_classic(base_size=11,base_family='Arial')+
  theme(plot.title=element_text(size=12,face='bold',hjust=.5),
        plot.subtitle=element_text(size=9,hjust=.5,margin=margin(b=6)),
        strip.background=element_rect(fill='#F5F5F5',colour='#777777',linewidth=.45),
        strip.text=element_text(size=11,face='italic'),
        panel.border=element_rect(fill=NA,colour='#777777',linewidth=.45),
        axis.line=element_blank(),axis.text=element_text(colour='#222222',size=10),
        axis.title.y=element_text(size=10),panel.spacing=grid::unit(.65,'lines'),
        legend.position='bottom',legend.title=element_text(size=8.5,face='bold'),
        legend.text=element_text(size=8.5),legend.key.size=grid::unit(.35,'cm'),
        legend.margin=margin(t=1),plot.margin=margin(7,7,4,7))
stem <- 'Figure6C_three_gene_GSM_boxplot_20260902'
ggsave(file.path(out,paste0(stem,'.png')),p,width=6.4,height=3.8,dpi=400,bg='white',device=grDevices::png,type='cairo')
ggsave(file.path(out,paste0(stem,'.pdf')),p,width=6.4,height=3.8,device=cairo_pdf,bg='white')
ggsave(file.path(out,paste0(stem,'.tif')),p,width=6.4,height=3.8,dpi=600,compression='lzw',bg='white',device=grDevices::tiff,type='cairo')
for(ext in c('png','pdf','tif')) {
  dest<-file.path(figure_dir,paste0(stem,'.',ext))
  if(file.exists(dest))stop('Refusing to overwrite existing figure: ',dest)
  stopifnot(file.copy(file.path(out,paste0(stem,'.',ext)),dest))
}

message('[4/4] Saving provenance, caption and checks')
caption<-paste(
  '(C) GSM-level pseudobulk expression of AGT, CXCR4, and FOS across all five revised fibroblast states.',
  'Raw RNA counts were summed within each GSM across qualified fibroblasts. Samples with at least 50 fibroblasts were retained (16 NC and 12 PD GSMs from GSE164241 and GSE171213).',
  'Pseudobulk libraries were normalized by the trimmed mean of M values (TMM), and expression is displayed as log2 counts per million with a prior count of 2.',
  'Each point represents one GSM, with color indicating dataset; NC and PD boxes show the median and interquartile range, with whiskers extending to observations within 1.5 times the interquartile range.',
  'Each gene has its own vertical scale. These unadjusted distributions are descriptive; no differential-expression significance is inferred from the box plots, and differences may reflect both within-state expression and fibroblast-state composition.')
writeLines(caption,file.path(out,'Figure6C_caption_en.txt'))
write.csv(data.frame(package=pkgs,version=vapply(pkgs,function(x)as.character(packageVersion(x)),character(1))),file.path(out,'package_versions.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
write.csv(data.frame(input=input,size=file.info(input)$size,modified=as.character(file.info(input)$mtime),md5=unname(tools::md5sum(input))),file.path(out,'input_provenance.csv'),row.names=FALSE)
writeLines(c(
  'PASS: source is the 16,032-cell, five-state Figure 6 fibroblast object.',
  'PASS: RNA counts are nonnegative integers; count columns match metadata rows.',
  'PASS: raw count totals are conserved in GSM aggregation.',
  'PASS: each GSM has one group/dataset label; 28 eligible GSMs (16 NC, 12 PD).',
  'PASS: exactly one plotted point per gene per eligible GSM (84 total points).',
  'No source Seurat object modified. No new R package installed. No cell-level significance tests.',
  'These distributions are not adjusted for dataset or fibroblast-state composition; use descriptive interpretation.'
),file.path(out,'verification.txt'))
print(sm);print(summ)
message('FIGURE6C COMPLETE')
