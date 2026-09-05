.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5","G:/gurobi/Rlib45",.libPaths()))
library(ggplot2)
fit<-readRDS("Locked_three_gene_logistic_model.rds")
genes<-c("AGT","CXCR4","FOS");beta<-coef(fit)[genes]
vals<-lapply(genes,function(g)seq(min(fit$model[[g]]),max(fit$model[[g]]),length.out=3));names(vals)<-genes
contrib<-lapply(genes,function(g)beta[g]*vals[[g]]);names(contrib)<-genes
mins<-vapply(contrib,min,numeric(1));ranges<-vapply(contrib,function(x)diff(range(x)),numeric(1));sc<-100/max(ranges)
mt<-sum(ranges)*sc;lp<-coef(fit)[1]+sum(mins)
rows<-do.call(rbind,lapply(seq_along(genes),function(i)data.frame(Gene=genes[i],Expression=vals[[i]],Points=(contrib[[i]]-mins[i])*sc,y=5-i)))
segs<-do.call(rbind,lapply(split(rows,rows$Gene),function(x)data.frame(x=min(x$Points),xe=max(x$Points),y=x$y[1])))
tt<-seq(0,mt,length.out=6);pr<-c(.1,.5,.9);pp<-(qlogis(pr)-lp)*sc;k<-pp>=0&pp<=mt
p<-ggplot(rows,aes(Points,y))+geom_segment(data=segs,aes(x=x,xend=xe,y=y,yend=y),inherit.aes=FALSE,linewidth=.6)+geom_point(color="#4056C9",size=1.7)+geom_text(aes(label=sprintf("%.1f",Expression)),vjust=-1,size=2.7)+
 annotate("segment",x=0,xend=mt,y=1.5,yend=1.5,linewidth=.7)+annotate("point",x=tt,y=1.5,size=1.4)+annotate("text",x=tt,y=1.5,label=sprintf("%.0f",tt),vjust=-1,size=2.7)+
 annotate("segment",x=0,xend=mt,y=.65,yend=.65,linewidth=.7)+annotate("segment",x=pp[k],xend=pp[k],y=.58,yend=.72)+annotate("text",x=pp[k],y=.40,label=sprintf("%.2f",pr[k]),size=2.7)+
 annotate("text",x=-.035*mt,y=c(4,3,2,1.5,.65),label=c(genes,"Total points","PD probability"),hjust=1,fontface="bold",size=3.2)+
 scale_x_continuous(position="top",breaks=tt,labels=sprintf("%.0f",tt),limits=c(-.35*mt,1.05*mt))+scale_y_continuous(limits=c(.1,4.6),breaks=NULL)+
 labs(title="Three-gene nomogram",subtitle="Coefficient-weighted expression points",x="Points",y=NULL)+theme_classic(base_size=12)+theme(plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,size=9),axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(color="black"),plot.margin=margin(10,10,10,10))
ggsave("panels/Figure5_G_nomogram.pdf",p,width=6.6,height=4.8,device=cairo_pdf,bg="white")
ggsave("panels/Figure5_G_nomogram.png",p,width=6.6,height=4.8,dpi=300,bg="white")
write.csv(rows,"nomogram_points_verified.csv",row.names=FALSE)
stopifnot(max(abs(plogis(lp+(rowSums(sapply(genes,function(g)(beta[g]*fit$model[[g]]-mins[g])*sc)))/sc)-fitted(fit)))<1e-8)
