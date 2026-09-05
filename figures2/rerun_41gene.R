options(stringsAsFactors=FALSE)
.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5", "G:/gurobi/Rlib45", .libPaths()))
pk <- c("caret","DALEX","ggplot2","randomForest","kernlab","pROC","gbm","nnet","glmnet","rpart")
stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(caret);library(DALEX);library(ggplot2);library(randomForest);library(kernlab);library(pROC)})
out <- "three_gene_revision_20260902"
dir.create(file.path(out,"models"),showWarnings=FALSE,recursive=TRUE)
old <- new.env(); load(".RData",envir=old)
genes <- scan(file.path(out,"candidate41.txt"),what="",quiet=TRUE)
stopifnot(length(genes)==41,all(genes %in% rownames(old$expression_matrix)))
d <- as.data.frame(t(old$expression_matrix[genes,,drop=FALSE]))
d$Type <- old$sample_classification[rownames(d),"title"]
stopifnot(!anyNA(d))
set.seed(123)
idx <- createDataPartition(y=d$Type,p=.7,list=FALSE)
stopifnot(identical(as.integer(idx),as.integer(old$inTrain)))
train_data <- d[idx,]; test_data <- d[-idx,]
write.csv(data.frame(Sample=rownames(d),Group=d$Type,Split=ifelse(seq_len(nrow(d)) %in% idx,"Training","Screening holdout")),file.path(out,"sample_split_41gene.csv"),row.names=FALSE)
control <- trainControl(method="repeatedcv",number=5,savePredictions=TRUE)
methods <- c(RF="rf",SVM="svmRadial",GLM="glm",GBM="gbm",KNN="knn",NNET="nnet",LASSO="glmnet",DT="rpart")
models <- list()
for(nm in names(methods)) {
 cat("Fitting",nm,"at",format(Sys.time()),"\n"); flush.console()
 args <- list(form=Type~.,data=train_data,method=unname(methods[nm]),trControl=control)
 if(nm=="SVM") args$prob.model <- TRUE
 if(nm=="GLM") args$family <- "binomial"
 if(nm=="NNET") args$trace <- FALSE
 if(nm=="GBM") args$verbose <- FALSE
 models[[nm]] <- do.call(caret::train,args)
 saveRDS(models[[nm]],file.path(out,"models",paste0(nm,"_41gene.rds")))
}
saveRDS(list(models=models,train=train_data,test=test_data,genes=genes,seed=123,rng_after_fit=.Random.seed),file.path(out,"eight_models_and_data.rds"))
pred <- do.call(rbind,lapply(names(models),function(nm) do.call(rbind,lapply(c("Training","Screening holdout"),function(s){x<-if(s=="Training")train_data else test_data; data.frame(Model=nm,Split=s,Sample=rownames(x),Group=x$Type,Y=as.integer(x$Type=="Treat"),Probability=predict(models[[nm]],x,type="prob")[["Treat"]])}))))
pred$AbsResidual <- abs(pred$Y-pred$Probability)
write.csv(pred,file.path(out,"eight_models_sample_predictions.csv"),row.names=FALSE)
metrics <- do.call(rbind,lapply(split(pred,interaction(pred$Model,pred$Split,drop=TRUE)),function(z)data.frame(Model=z$Model[1],Split=z$Split[1],N=nrow(z),AUC=as.numeric(pROC::auc(pROC::roc(z$Y,z$Probability,levels=c(0,1),direction="<",quiet=TRUE))),RMSE=sqrt(mean((z$Y-z$Probability)^2)))))
write.csv(metrics,file.path(out,"eight_models_metrics.csv"),row.names=FALSE)
print(metrics)
capture.output(lapply(models,function(m)m$bestTune),file=file.path(out,"best_tuning_parameters.txt"))
cat("Computing original DALEX permutation rankings\n")
importance <- list()
for(nm in c("RF","SVM","GLM","KNN","GBM","NNET","LASSO","DT")) {
 ex <- DALEX::explain(models[[nm]],label=nm,data=test_data,y=as.integer(test_data$Type=="Treat"),predict_function=function(object,newdata)predict(object,newdata,type="prob")[["Treat"]],verbose=FALSE)
 vi <- DALEX::variable_importance(ex)
 saveRDS(vi,file.path(out,"models",paste0(nm,"_permutation_importance.rds")))
 vv <- as.data.frame(vi)
 write.csv(vv,file.path(out,paste0(nm,"_importance_raw.csv")),row.names=FALSE)
 vv <- vv[vv$variable %in% genes,]
 ag <- aggregate(dropout_loss~variable,vv,mean)
 ag <- ag[order(-ag$dropout_loss,match(ag$variable,genes)),]
 ag$Rank <- seq_len(nrow(ag));ag$Model<-nm
 importance[[nm]]<-ag
}
im <- do.call(rbind,importance)
write.csv(im,file.path(out,"eight_models_importance_rankings.csv"),row.names=FALSE)
top <- lapply(importance,function(z)head(z$variable,20))
consensus <- Reduce(intersect,top[c("RF","SVM","NNET","LASSO")])
writeLines(consensus,file.path(out,"four_method_top20_consensus.txt"))
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
write.csv(data.frame(Package=pk,Version=vapply(pk,function(p)as.character(packageVersion(p)),character(1))),file.path(out,"R_package_versions.csv"),row.names=FALSE)
cat("DONE. Four-method top20 consensus:",paste(consensus,collapse=", "),"\n")
