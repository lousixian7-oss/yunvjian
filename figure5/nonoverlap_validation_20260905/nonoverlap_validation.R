try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE)
.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5','G:/gurobi/Rlib45',.libPaths()))
options(stringsAsFactors=FALSE)
for(pkg in c('affy','hgu133plus2cdf','hgu133plus2.db','AnnotationDbi','ggplot2')) stopifnot(requireNamespace(pkg,quietly=TRUE))
root <- normalizePath('tripod_20260905',winslash='/')
out <- file.path(root,'nonoverlap_validation')
archive <- 'G:/1Yunvjian/1A\u7389\u5973\u714e/code'
genes <- c('AGT','CXCR4','FOS')
meta <- read.csv(file.path(root,'overlap_audit/GSE16134_nonoverlapping_patient_samples.csv'))
stopifnot(nrow(meta)==66,length(unique(meta$patient))==30,sum(meta$group=='NC')==6)
meta$y <- as.integer(meta$group=='PD')
meta$patient_id <- paste0('GSE16134_P',meta$patient)
files <- file.path(out,'CEL',paste0(meta$sample_id,'.CEL.gz'))
stopifnot(all(file.exists(files)))
cat('Starting RMA on the 66 nonoverlapping arrays only.\n');flush.console()
if(!file.exists(file.path(out,'RMA66_probe_expression.rds'))){
 raw <- affy::ReadAffy(filenames=files)
 eset <- affy::rma(raw)
 probe <- Biobase::exprs(eset);colnames(probe)<-meta$sample_id
 saveRDS(probe,file.path(out,'RMA66_probe_expression.rds'))
 rm(raw,eset);gc()
} else probe<-readRDS(file.path(out,'RMA66_probe_expression.rds'))
stopifnot(identical(colnames(probe),meta$sample_id))
ann<-AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db,keys=rownames(probe),keytype='PROBEID',columns='SYMBOL')
ann<-ann[!is.na(ann$SYMBOL)&ann$SYMBOL!='',];ann<-ann[!duplicated(ann$PROBEID),]
selected<-ann[ann$SYMBOL %in% genes,]
write.csv(selected,file.path(out,'three_gene_probe_mapping.csv'),row.names=FALSE)
expr66<-sapply(genes,function(g) colMeans(probe[selected$PROBEID[selected$SYMBOL==g],,drop=FALSE]))
rownames(expr66)<-meta$sample_id
write.csv(expr66,file.path(out,'RMA66_three_gene_expression.csv'))
fit<-readRDS(file.path(archive,'coredata/model/Locked_three_gene_logistic_model.rds'))
expected<-c(-28.1615557900581,2.02555945620692,1.70618768640773,0.412233157693864)
stopifnot(max(abs(coef(fit)-expected))<1e-10)
dev<-t(as.matrix(read.csv(file.path(out,'development_three_genes.csv'),row.names=1,check.names=FALSE)))[,genes]
split<-read.csv(file.path(archive,'coredata/bulk/patient_split.csv'))
tr<-dev[split$sample_id[split$cohort=='Training'],,drop=FALSE]
stopifnot(nrow(tr)==120)
mu<-colMeans(tr);sdev<-apply(tr,2,sd);beta<-coef(fit)[genes];bz<-beta*sdev;iz<-unname(coef(fit)[1]+sum(beta*mu))
score<-function(m) plogis(iz+drop(scale(m[,genes,drop=FALSE])%*%bz))
primary<-score(expr66)
full<-t(as.matrix(read.csv(file.path(out,'full_cohort_three_genes.csv'),row.names=1,check.names=FALSE)))[,genes]
subset_rescale<-score(full[meta$sample_id,,drop=FALSE])
legacy<-read.csv(file.path(root,'probability_reconciliation/GSE16134_fixed_predictions_S6.csv'))
legacy_p<-legacy$probability[match(meta$sample_id,legacy$sample_id)]
stopifnot(max(abs(score(full)[match(meta$sample_id,rownames(full))]-legacy_p))<1e-10)
pred<-meta[,c('sample_id','patient_id','group','y')]
pred$RMA66_subset_scaling<-primary
pred$RMA310_subset_scaling<-subset_rescale
pred$RMA310_full_scaling_filtered<-legacy_p
stopifnot(all(is.finite(as.matrix(pred[,5:7]))))
write.csv(pred,file.path(out,'predictions_66.csv'),row.names=FALSE)
parameters<-data.frame(gene=genes,locked_beta=beta,training_mean=mu,training_sd=sdev,standardized_beta=bz,RMA66_mean=colMeans(expr66),RMA66_sd=apply(expr66,2,sd))
write.csv(parameters,file.path(out,'transport_parameters.csv'),row.names=FALSE)
writeLines(paste('Locked raw intercept:',coef(fit)[1],'\nEquivalent standardized intercept:',iz),file.path(out,'intercepts.txt'))
auc<-function(y,p){n1<-sum(y==1);n0<-sum(y==0);if(n1*n0==0)return(NA_real_);(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)}
metrics<-function(y,p,calibration=TRUE){
 lp<-qlogis(pmin(pmax(p,1e-8),1-1e-8))
 z<-c(AUC=auc(y,p),Brier=mean((p-y)^2),Sensitivity_0.5=mean(p[y==1]>=.5),Specificity_0.5=mean(p[y==0]<.5),Mean_prediction=mean(p))
 if(calibration){
 f1<-suppressWarnings(glm(y~1,offset=lp,family=binomial()));f2<-suppressWarnings(glm(y~lp,family=binomial()))
 z<-c(z,Calibration_intercept=if(f1$converged)unname(coef(f1)[1]) else NA_real_,Calibration_slope=if(f2$converged)unname(coef(f2)[2]) else NA_real_)
 }
 z
}
paths<-names(pred)[5:7]
point<-do.call(rbind,lapply(paths,function(k) data.frame(path=k,t(metrics(pred$y,pred[[k]])),check.names=FALSE)))
write.csv(point,file.path(out,'performance_comparison.csv'),row.names=FALSE)
print(point)
cat('Patient-cluster bootstrap: 2000 replicates, fixed predictions and preprocessing.\n');flush.console()
set.seed(20260905);B<-2000;patients<-unique(pred$patient_id)
boot<-matrix(NA_real_,B,7,dimnames=list(NULL,names(metrics(pred$y,primary))))
boot_comp<-matrix(NA_real_,B,4,dimnames=list(NULL,c('RMA310_subset_AUC','RMA310_subset_Brier','legacy_AUC','legacy_Brier')))
for(i in seq_len(B)){
 draw<-sample(patients,length(patients),replace=TRUE)
 idx<-unlist(lapply(draw,function(id)which(pred$patient_id==id)),use.names=FALSE)
 if(length(unique(pred$y[idx]))<2)next
 boot[i,]<-metrics(pred$y[idx],primary[idx])
 boot_comp[i,]<-c(auc(pred$y[idx],subset_rescale[idx]),mean((pred$y[idx]-subset_rescale[idx])^2),auc(pred$y[idx],legacy_p[idx]),mean((pred$y[idx]-legacy_p[idx])^2))
}
write.csv(boot,file.path(out,'patient_cluster_bootstrap_2000.csv'),row.names=FALSE)
write.csv(boot_comp,file.path(out,'sensitivity_bootstrap_2000.csv'),row.names=FALSE)
ci<-data.frame(metric=colnames(boot),estimate=as.numeric(metrics(pred$y,primary)),lower=apply(boot,2,function(v)quantile(v,.025,na.rm=TRUE)),upper=apply(boot,2,function(v)quantile(v,.975,na.rm=TRUE)),valid_replicates=colSums(is.finite(boot)))
write.csv(ci,file.path(out,'primary_metrics_95CI.csv'),row.names=FALSE)
lopo<-do.call(rbind,lapply(patients,function(id){k<-pred$patient_id!=id;data.frame(omitted_patient=id,AUC=auc(pred$y[k],primary[k]),Brier=mean((pred$y[k]-primary[k])^2))}))
write.csv(lopo,file.path(out,'leave_one_patient_out_influence.csv'),row.names=FALSE)
thresholds<-seq(.05,.95,.01);prev<-mean(pred$y)
dca<-data.frame(threshold=thresholds,model=sapply(thresholds,function(t)mean(primary>=t & pred$y==1)-mean(primary>=t & pred$y==0)*t/(1-t)),treat_all=prev-(1-prev)*thresholds/(1-thresholds),treat_none=0)
write.csv(dca,file.path(out,'decision_curve.csv'),row.names=FALSE)
# Descriptive four-bin calibration; no recalibration is applied to predictions.
bin<-cut(primary,unique(quantile(primary,seq(0,1,length.out=5))),include.lowest=TRUE)
cal<-aggregate(cbind(predicted=primary,observed=pred$y),list(bin=bin),mean)
cal$n<-as.integer(table(bin));write.csv(cal,file.path(out,'calibration_bins.csv'),row.names=FALSE)
library(ggplot2)
theme_set(theme_classic(base_size=12))
roc<-data.frame(FPR=c(0,sapply(sort(unique(primary),decreasing=TRUE),function(t)mean(primary[pred$y==0]>=t)),1),TPR=c(0,sapply(sort(unique(primary),decreasing=TRUE),function(t)mean(primary[pred$y==1]>=t)),1))
p1<-ggplot(roc,aes(FPR,TPR))+geom_step(color='#207EA4',linewidth=1)+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Participant-disjoint subset',subtitle=sprintf('66 sites / 30 participants; AUC %.3f (95%% CI %.3f-%.3f)',ci$estimate[1],ci$lower[1],ci$upper[1]),x='1 - Specificity',y='Sensitivity')
p2<-ggplot(cal,aes(predicted,observed))+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+geom_point(size=3,color='#207EA4')+geom_line(color='#207EA4')+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Descriptive calibration',subtitle='Four quantile bins; only six unaffected sites',x='Mean predicted probability',y='Observed affected-site proportion')
long<-rbind(data.frame(threshold=thresholds,NB=dca$model,strategy='Model'),data.frame(threshold=thresholds,NB=dca$treat_all,strategy='Treat all'),data.frame(threshold=thresholds,NB=0,strategy='Treat none'))
p3<-ggplot(long,aes(threshold,NB,color=strategy))+geom_line(linewidth=.8)+scale_color_manual(values=c('Model'='#207EA4','Treat all'='#D4844A','Treat none'='grey50'))+labs(title='Exploratory decision curve',subtitle='Net benefit reflects the selected tissue-sample composition',x='Threshold probability',y='Net benefit',color=NULL)+theme(legend.position='bottom')
plots<-list(ROC=p1,Calibration=p2,Decision_curve=p3)
for(n in names(plots)){ggsave(file.path(out,paste0(n,'.png')),plots[[n]],width=6.5,height=5,dpi=180,bg='white');ggsave(file.path(out,paste0(n,'.pdf')),plots[[n]],width=6.5,height=5)}
saveRDS(list(locked_coefficients=coef(fit),intercept_z=iz,parameters=parameters,bootstrap_seed=20260905,bootstrap_B=B),file.path(out,'analysis_parameters.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('FINISHED\n');print(ci)
