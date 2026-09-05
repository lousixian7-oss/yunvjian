options(stringsAsFactors=FALSE)
.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5','G:/gurobi/Rlib45',.libPaths()))
stopifnot(requireNamespace('ggplot2',quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2))
out<-'all_subset_model_audit_20260902/Figure6C_three_gene_GSM'
plotdata<-readRDS(file.path(out,'Figure6C_plot_data.rds'))
eligible<-read.csv(file.path(out,'Table_TMM_normalization.csv'))
eligible$group<-factor(eligible$group,levels=c('NC','PD'))
eligible$dataset<-factor(eligible$dataset)
plotdata$group<-factor(plotdata$group,levels=c('NC','PD'))
plotdata$dataset<-factor(plotdata$dataset,levels=levels(eligible$dataset))
plotdata$gene<-factor(plotdata$gene,levels=c('AGT','CXCR4','FOS'))
set.seed(20260902)
# Recreate the original plot directly from saved sample-level values, without
# loading the Seurat object or refitting any expression model.
code<-readLines('all_subset_model_audit_20260902/make_Figure6C_three_gene_GSM.R',warn=FALSE)
start<-grep('^group_cols<-',code);end<-grep('^stem <-',code)
stopifnot(length(start)==1L,length(end)==1L)
eval(parse(text=code[start:(end-1L)]))
stats<-read.csv(file.path(out,'Table_three_gene_dataset_adjusted_edgeR.csv'))
stats$gene<-factor(stats$gene,levels=levels(plotdata$gene))
expected<-ifelse(stats$BH_q_three_targets<.001,'***',ifelse(stats$BH_q_three_targets<.01,'**',ifelse(stats$BH_q_three_targets<.05,'*','ns')))
stopifnot(identical(expected,stats$significance),identical(stats$significance,c('**','ns','ns')))
ann<-do.call(rbind,lapply(seq_len(nrow(stats)),function(i){
  vals<-plotdata$log2_CPM[plotdata$gene==stats$gene[i]]
  span<-max(diff(range(vals)),1)
  data.frame(gene=stats$gene[i],y=max(vals)+.10*span,tip=.035*span,
             text_y=max(vals)+.19*span,label=stats$significance[i])
}))
p<-p+
  geom_segment(data=ann,aes(x=1,xend=2,y=y,yend=y),inherit.aes=FALSE,linewidth=.5,colour='#333333')+
  geom_segment(data=ann,aes(x=1,xend=1,y=y,yend=y-tip),inherit.aes=FALSE,linewidth=.5,colour='#333333')+
  geom_segment(data=ann,aes(x=2,xend=2,y=y,yend=y-tip),inherit.aes=FALSE,linewidth=.5,colour='#333333')+
  geom_text(data=ann,aes(x=1.5,y=text_y,label=label),inherit.aes=FALSE,size=4.5,fontface='bold')+
  labs(caption='Dataset-adjusted tests; BH correction across three genes')+
  theme(plot.caption=element_text(size=8,hjust=.5,margin=margin(t=2)))
stem<-'Figure6C_three_gene_GSM_boxplot_significance_20260902'
ggsave(file.path(out,paste0(stem,'.png')),p,width=6.4,height=4.0,dpi=400,bg='white',device=grDevices::png,type='cairo')
ggsave(file.path(out,paste0(stem,'.pdf')),p,width=6.4,height=4.0,device=cairo_pdf,bg='white')
ggsave(file.path(out,paste0(stem,'.tif')),p,width=6.4,height=4.0,dpi=600,bg='white',device=grDevices::tiff,type='cairo',compression='lzw')
figure_dir<-'G:/1Yunvjian/1A\u7389\u5973\u714e/Figure'
for(ext in c('png','pdf','tif'))stopifnot(file.copy(file.path(out,paste0(stem,'.',ext)),file.path(figure_dir,paste0(stem,'.',ext)),overwrite=TRUE))
write.csv(ann,file.path(out,'Figure6C_significance_annotations.csv'),row.names=FALSE)
cat('Saved annotation-only revision; no values or statistical results changed.\n')
