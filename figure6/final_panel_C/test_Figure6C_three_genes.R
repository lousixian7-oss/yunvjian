options(stringsAsFactors=FALSE)
.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5','G:/gurobi/Rlib45',.libPaths()))
stopifnot(requireNamespace('edgeR',quietly=TRUE))
suppressPackageStartupMessages(library(edgeR))
out<-'all_subset_model_audit_20260902/Figure6C_three_gene_GSM'
pb<-readRDS(file.path(out,'GSM_pseudobulk_counts_and_metadata.rds'))
sm<-pb$sample_info
sm$group<-factor(sm$group,levels=c('NC','PD'))
sm$dataset<-droplevels(factor(sm$dataset))
design<-model.matrix(~dataset+group,data=sm)
stopifnot(qr(design)$rank==ncol(design),identical(colnames(pb$counts),sm$GSM))
norm<-read.csv(file.path(out,'Table_TMM_normalization.csv'))
stopifnot(identical(norm$GSM,sm$GSM))
y<-DGEList(counts=pb$counts,lib.size=norm$library_size,norm.factors=norm$TMM_factor)
auto_keep<-filterByExpr(y,design=design)
target_index<-match(pb$targets,rownames(y));stopifnot(!anyNA(target_index))
# Retain the three already-plotted targets for the requested targeted tests.
# Report if any would have failed the usual expression filter.
keep<-auto_keep;keep[target_index]<-TRUE
y<-y[keep,,keep.lib.sizes=TRUE]
fit<-glmQLFit(y,design,robust=TRUE)
test<-glmQLFTest(fit,coef='groupPD')
full<-topTags(test,n=Inf,sort.by='none')$table
res<-data.frame(gene=pb$targets,log2FC_PD_vs_NC=full[pb$targets,'logFC'],
               logCPM=full[pb$targets,'logCPM'],QL_F=full[pb$targets,'F'],
               P_value=full[pb$targets,'PValue'],
               FDR_all_tested_genes=full[pb$targets,'FDR'],
               passes_standard_expression_filter=auto_keep[target_index],
               n_NC=sum(sm$group=='NC'),n_PD=sum(sm$group=='PD'),
               model='~ dataset + group',genes_tested=nrow(full))
res$BH_q_three_targets<-p.adjust(res$P_value,'BH')
res$significance<-ifelse(res$BH_q_three_targets<.001,'***',ifelse(res$BH_q_three_targets<.01,'**',ifelse(res$BH_q_three_targets<.05,'*','ns')))
res$direction<-ifelse(res$log2FC_PD_vs_NC>0,'PD higher','PD lower')
write.csv(res,file.path(out,'Table_three_gene_dataset_adjusted_edgeR.csv'),row.names=FALSE)
write.csv(data.frame(gene=rownames(full),full),file.path(out,'Table_all_tested_genes_edgeR.csv'),row.names=FALSE)
write.csv(data.frame(GSM=sm$GSM,design),file.path(out,'Table_edgeR_design_matrix.csv'),row.names=FALSE)
write.csv(as.data.frame(table(sm$dataset,sm$group)),file.path(out,'Table_test_samples_by_dataset_group.csv'),row.names=FALSE)
saveRDS(list(fit=fit,test=test,target_results=res),file.path(out,'Three_gene_edgeR_targeted_fit.rds'))
writeLines(c('Three targeted exploratory tests of AGT, CXCR4 and FOS in all revised fibroblasts.',
 'Source: same 28 GSMs and pseudobulk counts used in Figure 6C; 16 NC, 12 PD.',
 'Model: robust edgeR quasi-likelihood with dataset and group; PD vs NC.',
 'Library sizes and TMM factors retained exactly from the existing plot; no renormalization after filtering.',
 'filterByExpr used for dispersion estimation; preselected targets retained, with filter status reported.',
 'Primary correction: Benjamini-Hochberg across three requested targets. Transcriptome-wide FDR also reported.',
 'No pooling of individual cells as replicates. No change of cell population based on significance.',
 'These associations combine within-state expression changes and fibroblast-state compositional differences.',
 'A non-significant result does not establish absence of association.'),file.path(out,'Three_gene_tests_methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'Three_gene_tests_sessionInfo.txt'))
print(res,row.names=FALSE)
cat('TARGETED TESTS COMPLETE\n')
