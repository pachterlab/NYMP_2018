
library(scde)

#args <- commandArgs(TRUE)
#base_dir <- args[[1]]
base_dir <- '/home/lynnyi/single/6.3_sims_exp/'
write_dir <- file.path(base_dir, 'de')
tximp <- readRDS(file.path(base_dir, 'tximport.rds'))
counts <- tximp$counts

sg <- factor(c(rep(1, 105), rep(2, 105)))
names(sg) <- as.character(1:210)
#need to recast



counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
colnames(counts) <- as.character(1:210)

runSCDE <- function(counts, sg) {
	error_models <- scde.error.models(counts = counts, groups = sg, n.cores = 20, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
	#saveRDS(error_models, file.path(write_dir, 'scde_error_models.rds'))
	prior <- scde.expression.prior(models=error_models, counts, show.plot=FALSE)
	scde_DE <- scde.expression.difference(error_models, counts, prior, sg, verbose=1, n.cores=20, n.randomizations=100)
	scde_DE$pvalue <- 2*pnorm(-abs(scde_DE$Z))
	#saveRDS(scde_DE, file.path(write_dir, 'scde.rds'))
}


time <- system.time(runSCDE(counts, sg))
saveRDS(time, '/home/lynnyi/single/R/scdetime.R')

