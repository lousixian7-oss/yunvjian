.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5', 'G:/gurobi/Rlib45', .libPaths()))
stopifnot(requireNamespace('CellChat', quietly=TRUE), requireNamespace('ggplot2', quietly=TRUE))
suppressPackageStartupMessages(library(CellChat))
library(ggplot2)
objects <- lapply(c('NC','PD'), function(g) readRDS(paste0('inputs/',g,'_cellchat.rds')))
names(objects) <- c('NC','PD')
sink('saved_model_parameters.txt')
for (g in names(objects)) {
  cat('\n',g,'\n'); print(objects[[g]]@options); print(table(objects[[g]]@idents))
}
cat('\nDatabase identical:', identical(objects$NC@DB, objects$PD@DB), '\n')
sink()
targets <- c('Reticular stroma','Immune stroma','Basal epithelium','Suprabasal epithelial','Inflamed epithelium')
# Extract directly from saved arrays, including zero / nonsignificant results.
tables <- lapply(names(objects), function(g) {
  obj <- objects[[g]]
  lr <- obj@LR$LRsig
  lr <- lr[lr$pathway_name == 'TGFb', , drop=FALSE]
  rows <- lapply(intersect(rownames(lr), dimnames(obj@net$prob)[[3]]), function(k) {
    data.frame(dataset=g, source='Reticular stroma', target=targets,
      interaction_name=k, interaction_name_2=lr[k,'interaction_name_2'],
      ligand=lr[k,'ligand'], receptor=lr[k,'receptor'],
      prob=as.numeric(obj@net$prob['Reticular stroma',targets,k]),
      pval=as.numeric(obj@net$pval['Reticular stroma',targets,k]))
  })
  do.call(rbind,rows)
})
all_results <- do.call(rbind,tables)
all_results$passes_nominal_p <- all_results$prob > 0 & all_results$pval < 0.05
write.csv(all_results, 'TGFb_ReticularStroma_all_tested.csv',row.names=FALSE)
df <- all_results[all_results$passes_nominal_p, ]
write.csv(df, 'TGFb_ReticularStroma_plotted.csv',row.names=FALSE)
print(table(df$dataset)); print(df)
stopifnot(nrow(df)>0, all(df$pval<0.05), all(df$prob>0))
shown_targets <- targets[targets %in% df$target]
df$target <- factor(df$target,levels=shown_targets)
df$dataset <- factor(df$dataset,levels=c('NC','PD'))
pairs <- sort(unique(df$interaction_name_2))
df$pair <- factor(df$interaction_name_2, levels=rev(pairs))
df$p_category <- factor(ifelse(df$pval<=0.01,'P <= 0.01','0.01 < P < 0.05'),
  levels=c('0.01 < P < 0.05','P <= 0.01'))
target_labels <- c('Reticular stroma'='Reticular\nstroma','Immune stroma'='Immune\nstroma',
 'Basal epithelium'='Basal\nepithelium','Suprabasal epithelial'='Suprabasal\nepithelial','Inflamed epithelium'='Inflamed\nepithelium')
p <- ggplot(df,aes(x=dataset,y=pair)) +
 geom_point(aes(colour=prob,size=p_category)) +
 facet_grid(~target,drop=FALSE,labeller=as_labeller(target_labels)) +
 scale_x_discrete(drop=FALSE,expand=expansion(add=0.6)) +
 scale_y_discrete(expand=expansion(add=0.65)) +
 scale_colour_gradientn(colours=c('#5E4FA2','#3288BD','#66C2A5','#FFFFBF','#FDAE61','#D53E4F','#9E0142'),
   transform='log10',breaks=c(1e-4,3e-4,1e-3,3e-3), labels=c('0.0001','0.0003','0.001','0.003'),
   name='Communication\nprobability\n(log scale)') +
 scale_size_manual(values=c(2.7,4.3),drop=FALSE,name='Nominal P') +
 labs(title='TGF\u03b2 ligand-receptor interactions',
   subtitle='Sender: Reticular stroma | Receivers shown above',x=NULL,y=NULL,
   caption='Dots: within-condition CellChat P < 0.05; blank: criterion not met / not inferred.\nNo significant Reticular stroma self-interaction. P values are not NC-PD differential tests.') +
 theme_bw(base_size=12,base_family='Arial') +
 theme(panel.grid.major=element_line(colour='#EEEEEE',linewidth=0.35),panel.grid.minor=element_blank(),
   panel.border=element_rect(colour='#999999',linewidth=0.45),
   strip.background=element_rect(fill='#F5F6F7',colour='#CCCCCC',linewidth=0.4),
   strip.text=element_text(size=11,margin=margin(5,2,5,2)),
   axis.text=element_text(colour='#222222'),axis.text.y=element_text(size=11),
   axis.text.x=element_text(size=12,face='bold'),axis.ticks=element_blank(),
   panel.spacing=grid::unit(0.12,'cm'),legend.title=element_text(size=10.5),
   legend.text=element_text(size=10),legend.key.height=grid::unit(0.45,'cm'),
   plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=11),
   plot.caption=element_text(size=8.6,hjust=0,colour='#555555',margin=margin(t=9)),
   plot.title.position='plot',plot.caption.position='plot',plot.margin=margin(8,10,8,8)) +
 guides(colour=guide_colourbar(order=1,barheight=grid::unit(2.6,'cm')),size=guide_legend(order=2))
ggsave('Figure9L_TGFb_LR_NC_PD.pdf',p,device=cairo_pdf,width=8.9,height=4.45,units='in',bg='white')
saveRDS(p,'Figure9L_TGFb_LR_plot.rds')
writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
write.csv(data.frame(package=c('CellChat','ggplot2'), version=vapply(c('CellChat','ggplot2'),function(x) as.character(packageVersion(x)), character(1))), 'package_versions.csv',row.names=FALSE)
