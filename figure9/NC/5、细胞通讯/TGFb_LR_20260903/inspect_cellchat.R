.libPaths(c('C:/Users/32266/AppData/Local/R/win-library/4.5', 'G:/gurobi/Rlib45', .libPaths()))
stopifnot(requireNamespace('CellChat', quietly=TRUE))
library(CellChat)
cat('CellChat version:', as.character(packageVersion('CellChat')), '\n')
print(args(netVisual_bubble))
for (grp in c('NC', 'PD')) {
  e <- new.env()
  nm <- load(paste0('inputs/', grp, '_cellchat.rda'), envir=e)
  cat('\nGROUP', grp, 'OBJECTS', nm, '\n')
  objs <- mget(nm, envir=e)
  obj <- objs[[which(vapply(objs, inherits, logical(1), 'CellChat'))[1]]]
  print(levels(obj@idents))
  print(names(obj@net))
  dat <- subsetCommunication(obj, signaling='TGFb', thresh=1)
  dat$dataset <- grp
  write.csv(dat, paste0(grp, '_TGFb_all_existing_results.csv'), row.names=FALSE)
  print(dat[, intersect(c('source','target','ligand','receptor','prob','pval','interaction_name_2'), names(dat))])
  saveRDS(obj, paste0('inputs/', grp, '_cellchat.rds'))
}
writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
