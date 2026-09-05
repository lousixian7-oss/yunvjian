.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5',.libPaths()))
stopifnot(requireNamespace('SeuratObject',quietly=TRUE))
out<-'.code_release/yunvjian/coredata/metadata';dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(nm in c('global','fibroblast')) {
 src<-if(nm=='global') 'G:/1Yunvjian/0a26.7.7singlecell/注释26.8.6/yjsl_clean_allGSE.rds' else 'G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/results/objects/fibroblast_revision_final_no_unassigned.rds'
 cat('Reading',nm,'\n');flush.console();o<-readRDS(src)
 md<-o@meta.data;md$cell_id<-rownames(md)
 write.csv(md,gzfile(file.path(out,paste0(nm,'_full_metadata.csv.gz'))),row.names=FALSE)
 for(red in names(o@reductions)) {
  emb<-o@reductions[[red]]@cell.embeddings
  if(ncol(emb)>50)emb<-emb[,1:50,drop=FALSE]
  write.csv(data.frame(cell_id=rownames(emb),emb,check.names=FALSE),gzfile(file.path(out,paste0(nm,'_',red,'_embeddings.csv.gz'))),row.names=FALSE)
 }
 writeLines(c(paste('source:',src),paste('cells:',nrow(md)),paste('reductions:',paste(names(o@reductions),collapse=', '))),file.path(out,paste0(nm,'_object_source.txt')))
 cat(nm,nrow(md),'exported\n');rm(o,md);gc()
}
capture.output(sessionInfo(),file=file.path(out,'export_sessionInfo.txt'))
