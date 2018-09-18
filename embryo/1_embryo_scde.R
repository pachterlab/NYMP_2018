
# Perform SCDE on embryo dataset -----------------------------
library(scde)

readRDS('~/NYMP_2018/embryo/embryonic_gene_tpms.rds') -> genetpms
readRDS('~/NYMP_2018/embryo/embryonic_metadata.rds') -> coldata
sg <- factor(coldata$day)
names(sg) <- as.character(coldata$sample)

genetpms<-apply(genetpms,2,function(x) {storage.mode(x) <- 'integer'; x})
colnames(genetpms) <- names(sg)

error_models <- scde.error.models(counts = genetpms, groups = sg, n.cores = 20, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
saveRDS(error_models, file.path(write_dir, 'scde_error_models.rds'))
prior <- scde.expression.prior(models=error_models, genetpms, show.plot=FALSE)
scde_DE <- scde.expression.difference(error_models, genetpms, prior, sg, verbose=1, n.cores=20, n.randomizations=100)
scde_DE$pvalue <- 2*pnorm(-abs(scde_DE$Z))
saveRDS(scde_DE, file.path(write_dir, '~/NYMP_2018/embryo/scde.rds'))


