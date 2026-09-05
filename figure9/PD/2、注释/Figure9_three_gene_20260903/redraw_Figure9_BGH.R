options(stringsAsFactors=FALSE)
.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5","G:/gurobi/Rlib45",.libPaths()))
pk<-c("Seurat","SeuratObject","ggplot2","patchwork","RColorBrewer","Matrix")
stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(patchwork)})
dir.create("panels",showWarnings=FALSE);dir.create("single_gene",showWarnings=FALSE)
objects<-list()
for(gr in c("PD","NC")) {
 e<-new.env();load(file.path("inputs",paste0(gr,"_annotation.rda")),e)
 o<-e$stRNA;DefaultAssay(o)<-"Spatial"
 # Match the original annotation plotting script. Work in memory only.
 o<-NormalizeData(o,assay="Spatial",normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
 objects[[gr]]<-o
}
model_genes<-c("AGT","CXCR4","FOS")
markers<-c("PI16","COL1A1","FN1","POSTN")
map_genes<-c(model_genes,markers)
bubble_genes<-c(model_genes,markers,"CD44")
stopifnot(all(vapply(objects,function(o)all(bubble_genes %in% rownames(o[["Spatial"]])),logical(1))))
expr<-lapply(objects,function(o)LayerData(o,assay="Spatial",layer="data")[bubble_genes,,drop=FALSE])
limits<-setNames(vapply(map_genes,function(g){x<-unlist(lapply(expr,function(z)as.numeric(z[g,])));pos<-x[x>0];if(length(pos))as.numeric(quantile(pos,.95))else 1},numeric(1)),map_genes)
write.csv(data.frame(Gene=map_genes,Min=0,Max=limits,Rule="Pooled PD+NC positive-expression 95th percentile; shared for the same gene"),"spatial_color_limits.csv",row.names=FALSE)
palette<-rev(RColorBrewer::brewer.pal(11,"Spectral"))
for(gr in names(objects)) {
 plots<-lapply(map_genes,function(g){
  pp<-SpatialFeaturePlot(objects[[gr]],features=g,pt.size.factor=5,min.cutoff=0,max.cutoff=limits[g],combine=FALSE)[[1]]+
   scale_fill_gradientn(colors=palette,limits=c(0,limits[g]),breaks=c(0,limits[g]/2,limits[g]),labels=function(x)sprintf("%.1f",x),oob=scales::squish,name=g)+
   guides(fill=guide_colorbar(direction="horizontal",title.position="top",title.hjust=.5,barwidth=grid::unit(2.0,"cm"),barheight=grid::unit(.22,"cm")))+
   theme(legend.position="top",legend.title=element_text(size=11,face="bold"),legend.text=element_text(size=7),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),plot.margin=margin(1,1,1,1))
  ggsave(file.path("single_gene",paste0(gr,"_",g,".pdf")),pp,width=2.35,height=3.55,device=cairo_pdf,bg="white")
  pp
 })
 panel<-wrap_plots(plots,nrow=1)
 stem<-if(gr=="PD")"Figure9B_PD_three_gene_with_markers"else"Figure9G_NC_three_gene_with_markers"
 ggsave(file.path("panels",paste0(stem,".pdf")),panel,width=15.8,height=3.7,device=cairo_pdf,bg="white")
 # PNG/TIFF are rendered from the final PDF by export_and_check.py.
}
bubble<-do.call(rbind,lapply(names(expr),function(gr){z<-expr[[gr]];data.frame(Group=gr,Gene=rownames(z),Average_expression=Matrix::rowMeans(z),Percent_expressing=Matrix::rowMeans(z>0)*100,N_spots=ncol(z))}))
write.csv(bubble,"Figure9H_spot_level_expression_summary.csv",row.names=FALSE)
bubble$Gene<-factor(bubble$Gene,levels=bubble_genes);bubble$Group<-factor(bubble$Group,levels=c("NC","PD"))
ph<-ggplot(bubble,aes(Gene,Group))+geom_point(aes(size=Percent_expressing,color=Average_expression))+
 scale_color_viridis_c(option="plasma",name="Mean expression\n(log-normalized)")+
 scale_size_area(max_size=10,limits=c(0,100),breaks=c(25,50,75,100),name="Spots expressing (%)")+
 scale_y_discrete(expand=expansion(add=.65))+labs(x=NULL,y=NULL)+theme_classic(base_size=13)+
 theme(axis.text=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.box="horizontal",legend.title=element_text(size=9),legend.text=element_text(size=8),plot.margin=margin(10,10,10,10))+
 guides(color=guide_colorbar(direction="horizontal",title.position="top",barwidth=grid::unit(3,"cm"),barheight=grid::unit(.3,"cm")),size=guide_legend(nrow=1,title.position="top",label.position="bottom"))
ggsave("panels/Figure9H_PD_NC_three_gene_with_markers.pdf",ph,width=6.6,height=3.7,device=cairo_pdf,bg="white")
write.csv(data.frame(Group=names(objects),Spots=vapply(objects,ncol,numeric(1)),Source=c("PD annotation object","NC annotation object")),"annotation_object_spot_counts.csv",row.names=FALSE)
capture.output(sessionInfo(),file="R_sessionInfo.txt")
write.csv(data.frame(Package=pk,Version=vapply(pk,function(p)as.character(packageVersion(p)),character(1))),"R_package_versions.csv",row.names=FALSE)
cat("DONE: three replacement panels; raw annotation objects not modified.\n")
