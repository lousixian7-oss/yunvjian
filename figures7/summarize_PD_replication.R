options(width=180)
.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5','G:/gurobi/Rlib45',.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root <- 'G:/1a空转rperio'
samples <- c('GSM6258255','GSM6258256','GSM6258257','GSM6258258')
dirs <- c('G:/1a空转r/3、反卷机(RCTD)',file.path(root,'3、反卷机(RCTD) - 副本'),file.path(root,'GSM6258257'),file.path(root,'GSM6258258'))
set.seed(20260824)
allstats <- list();audits <- list();plots <- list();blocks <- list()
for(i in seq_along(samples)) {
  fs <- list.files(dirs[i],pattern='weights_full[.]csv$',full.names=TRUE)
  if(length(fs)!=1) { warning('Cannot uniquely identify weights for ',samples[i]);next }
  w <- read.csv(fs,check.names=FALSE)
  stopifnot(all(c('spot_id','ECM_Fib','Activated_Fib') %in% names(w)))
  w <- w[rowSums(w[,-1,drop=FALSE])>0,,drop=FALSE]
  pp <- cor.test(w$ECM_Fib,w$Activated_Fib,method='pearson')
  sp <- suppressWarnings(cor.test(w$ECM_Fib,w$Activated_Fib,method='spearman',exact=FALSE))
  allstats[[i]] <- data.frame(sample=samples[i],group=if(i==1)'NC' else 'PD',n_spots=nrow(w),pearson_r=unname(pp$estimate),pearson_p_spot_naive=pp$p.value,spearman_rho=unname(sp$estimate),spearman_p_spot_naive=sp$p.value)
  p <- ggplot(w,aes(ECM_Fib,Activated_Fib))+geom_point(size=1.8,alpha=.65,color='#2878B5')+geom_smooth(method='lm',formula=y~x,se=TRUE,color='#C82423')+theme_classic(base_size=13)+theme(aspect.ratio=1)+labs(title=paste(samples[i],if(i==1)'NC' else 'PD'),subtitle=sprintf('n = %d spots; r = %.3f; rho = %.3f',nrow(w),pp$estimate,sp$estimate),x='ECM_Fib fraction',y='Activated_Fib fraction')
  plots[[i]] <- p
  if(i>=3) {
    standalone <- p+labs(title=paste(samples[i],': ECM_Fib and Activated_Fib'),subtitle=sprintf('n = %d spots; Pearson r = %.3f (P = %.3g)\nSpearman rho = %.3f (P = %.3g)',nrow(w),pp$estimate,pp$p.value,sp$estimate,sp$p.value),caption='Continuous RCTD weights_full. P values are not adjusted for spatial autocorrelation.')+theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=11),plot.caption=element_text(size=9))
    ggsave(file.path(dirs[i],'PD_ECM_Fib_Activated_Fib_spot_correlation.pdf'),standalone,width=7,height=6)
  }
  # Fixed spatial-block bootstrap: uncertainty sensitivity, not a donor-level test.
  candidates <- list.files(dirs[i],pattern='^stRNA_RCTD_refinedFib_ECM.*[.]rds$',full.names=TRUE)
  target <- if(i==1) candidates[grepl('ECM[.]rds$',candidates)] else candidates[grepl('ECM_PD[.]rds$',candidates)]
  if(length(target)!=1) next
  obj <- readRDS(target)
  ca <- attributes(attributes(attr(obj,'images')[[1]])$boundaries[['centroids']])
  xy <- as.matrix(ca$coords);rownames(xy)<-ca$cells;xy<-xy[w$spot_id,,drop=FALSE]
  rm(obj);gc()
  dd<-as.matrix(dist(xy));diag(dd)<-Inf;nn<-median(apply(dd,1,min))
  audits[[i]]<-data.frame(sample=samples[i],median_nearest_distance=nn,mean_neighbors50=mean(rowSums(dd<50)),isolated_fraction50=mean(rowSums(dd<50)==0))
  for(mult in c(2,4)) {
    block <- interaction(floor((xy[,1]-min(xy[,1]))/(nn*mult)),floor((xy[,2]-min(xy[,2]))/(nn*mult)),drop=TRUE)
    ix <- split(seq_len(nrow(w)),block)
    boot <- replicate(1000,{jj<-unlist(ix[sample(seq_along(ix),length(ix),replace=TRUE)],use.names=FALSE);c(Pearson=cor(w$ECM_Fib[jj],w$Activated_Fib[jj]),Spearman=cor(w$ECM_Fib[jj],w$Activated_Fib[jj],method='spearman'))})
    ci <- t(apply(boot,1,quantile,probs=c(.025,.975),na.rm=TRUE))
    blocks[[paste(i,mult)]]<-data.frame(sample=samples[i],block_width_nn=mult,n_blocks=length(ix),metric=rownames(ci),lower=ci[,1],upper=ci[,2],B=1000,note='Spatial-block resampling sensitivity; no independent-donor inference')
  }
}
tab<-do.call(rbind,allstats);bt<-do.call(rbind,blocks)
bt$estimate <- ifelse(bt$metric=='Pearson',tab$pearson_r[match(bt$sample,tab$sample)],tab$spearman_rho[match(bt$sample,tab$sample)])
bp <- ggplot(bt,aes(estimate,sample,color=factor(block_width_nn)))+geom_vline(xintercept=0,linetype=2,color='grey50')+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',position=position_dodge(width=.5),width=.18)+geom_point(position=position_dodge(width=.5),size=2.5)+facet_wrap(~metric)+theme_classic(base_size=13)+labs(title='Spatial-block bootstrap sensitivity',subtitle='1,000 resamples; 95% percentile intervals',x='Correlation coefficient',y=NULL,color='Block width\n(spot spacing)',caption='Exploratory within-section uncertainty; not a patient-level comparison. Fixed non-overlapping grids.')
combined <- wrap_plots(plots,ncol=2)+plot_annotation(title='Within-dataset spatial replication',subtitle='Same refined reference and continuous RCTD weights_full',caption='Each panel is one section. Spot-level correlations are descriptive; donor independence is unverified.')
for(sample in samples[3:4]) {
  od<-file.path(root,sample)
  write.csv(tab,file.path(od,'cross_section_correlation_summary.csv'),row.names=FALSE)
  write.csv(bt,file.path(od,'spatial_block_bootstrap_sensitivity.csv'),row.names=FALSE)
  write.csv(do.call(rbind,audits),file.path(od,'cross_section_coordinate_audit.csv'),row.names=FALSE)
  ggsave(file.path(od,'PD_replication_comparison.pdf'),combined,width=11,height=11)
  ggsave(file.path(od,'PD_replication_comparison.png'),combined,width=11,height=11,dpi=180,bg='white')
  ggsave(file.path(od,'spatial_block_bootstrap_sensitivity.pdf'),bp,width=10,height=5.5)
}
print(tab);print(bt)
