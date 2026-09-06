try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE)
.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5","G:/gurobi/Rlib45",.libPaths()))
out <- normalizePath("tripod_20260905/manuscript_nonoverlap_update/updated",winslash="/")
data_dir <- "tripod_20260905/nonoverlap_validation"
pred<-read.csv(file.path(data_dir,"predictions_66.csv"));primary<-pred$RMA66_subset_scaling
ci<-read.csv(file.path(data_dir,"primary_metrics_95CI.csv"));cal<-read.csv(file.path(data_dir,"calibration_bins.csv"));dca<-read.csv(file.path(data_dir,"decision_curve.csv"));thresholds<-dca$threshold
library(ggplot2)
theme_set(theme_classic(base_size=12))
roc<-data.frame(FPR=c(0,sapply(sort(unique(primary),decreasing=TRUE),function(t)mean(primary[pred$y==0]>=t)),1),TPR=c(0,sapply(sort(unique(primary),decreasing=TRUE),function(t)mean(primary[pred$y==1]>=t)),1))
p1<-ggplot(roc,aes(FPR,TPR))+geom_step(color='#207EA4',linewidth=1)+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Participant-disjoint subset',subtitle=sprintf('66 sites / 30 participants; AUC %.3f (95%% CI %.3f-%.3f)',ci$estimate[1],ci$lower[1],ci$upper[1]),x='1 - Specificity',y='Sensitivity')
p2<-ggplot(cal,aes(predicted,observed))+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+geom_point(size=3,color='#207EA4')+geom_line(color='#207EA4')+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Descriptive calibration',subtitle='Four quantile bins; only six unaffected sites',x='Mean predicted probability',y='Observed affected-site proportion')
long<-rbind(data.frame(threshold=thresholds,NB=dca$model,strategy='Model'),data.frame(threshold=thresholds,NB=dca$treat_all,strategy='Treat all'),data.frame(threshold=thresholds,NB=0,strategy='Treat none'))
p3<-ggplot(long,aes(threshold,NB,color=strategy))+geom_line(linewidth=.8)+scale_color_manual(values=c('Model'='#207EA4','Treat all'='#D4844A','Treat none'='grey50'))+labs(title='Exploratory decision curve',subtitle='Net benefit reflects the selected tissue-sample composition',x='Threshold probability',y='Net benefit',color=NULL)+theme(legend.position='bottom')
plots<-list(ROC=p1,Calibration=p2,Decision_curve=p3)

p1<-p1+labs(title="ROC curve",subtitle=sprintf("AUC %.3f (95%% CI %.3f-%.3f)",ci$estimate[1],ci$lower[1],ci$upper[1]))
p2<-p2+labs(subtitle="Brier 0.205; intercept 4.011; slope 1.727")
p3<-p3+labs(subtitle="Selected tissue-sample composition")
draw<-function(){
 grid::grid.newpage()
 grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,3,heights=grid::unit(c(.45,4.65),"in"))))
 grid::grid.text("GSE16134 participant-disjoint subset: 66 sites from 30 patients (60 affected and 6 unaffected)",vp=grid::viewport(layout.pos.row=1,layout.pos.col=1:3),gp=grid::gpar(fontsize=13,fontface="bold"))
 for(i in 1:3){
  z<-list(p1,p2,p3)[[i]]+theme(plot.title=element_text(size=12),plot.subtitle=element_text(size=9),plot.margin=margin(9,9,9,16))
  print(z,vp=grid::viewport(layout.pos.row=2,layout.pos.col=i))
  grid::grid.text(LETTERS[i],x=grid::unit((i-1)/3+.01,"npc"),y=grid::unit(.865,"npc"),gp=grid::gpar(fontsize=15,fontface="bold"))
 }
 grid::popViewport()
}
png(file.path(out,"FigureS10.png"),width=16,height=5.1,units="in",res=300,type="cairo");draw();dev.off()
pdf(file.path(out,"FigureS10.pdf"),width=16,height=5.1);draw();dev.off()
tiff(file.path(out,"FigureS10.tif"),width=16,height=5.1,units="in",res=300,compression="lzw",type="cairo");draw();dev.off()
cat("Saved Figure S10 PNG/PDF/TIF.\n")
