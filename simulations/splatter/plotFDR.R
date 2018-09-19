
# Supp Fig 5a
source('~/single/R/roc_curve.R')

mast <- readRDS('./mast_log.rds')
deseq2 <- readRDS('./deseq2.rds')
scde <- readRDS('./scde.rds')

#seurat <- readRDS('./seurat.rds')
tobit <- readRDS('./tobit.rds')
LR_univariate <- readRDS('./LR_univariate_deseq.rds')
LR_tx <- readRDS('./LR_deseq.rds')

mast$mast <- mast$p_val
deseq2$deseq2 <- deseq2$p_val
tobit$tobit <- tobit$p_val
scde$scde <- scde$pvalue
LR_univariate$LR_univariate <- LR_univariate$pvals
LR_tx$LR_tx <- LR_tx$pvalues

columns <- c('deseq2', 'tobit', 'mast', 'scde', 'LR_tx', 'LR_univariate')
names <- c('DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE', 'log reg - transcripts', 'log reg - gene counts')

results <- list(deseq2, tobit, mast, scde)
results <- lapply(results, function(x) { x$genes <- rownames(x)
					return(x) })
results <- append(results, list(LR_tx, LR_univariate))
table <- merge(results[[1]], results[[2]], by='genes', all=TRUE)
for(i in 3:length(results))
{
	table <- merge(table, results[[i]], by='genes', all=TRUE)
}
perturb <- read.table('./perturb.tsv', header=TRUE, stringsAsFactors=FALSE)
t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
colnames(t2g) <- c('transcripts', 'genes')

library(dplyr)
perturb <- merge(perturb, t2g, by.x= 'perturbed_transcripts', by.y = 'transcripts')
perturb <- perturb %>% group_by(genes) %>% summarise(perturb= any(perturbed))
table <- merge(table, perturb, by='genes', all=TRUE)
fdrs <- lapply(columns, function(x) calculate_fdr(table$perturb, table[[x]], x))


fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sensitivity, FDR = fdrs[[x]]$fdr, Method = names[[x]]))
fdrs <- do.call(rbind, fdrs)

library(ggplot2)
tiff('./fdr.tiff', height= 3000, width= 4000, res=500)
g <- ggplot(fdrs, aes(FDR, sensitivity, colour=Method)) + geom_path() + scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue', 'purple', 'purple4')) + theme_grey(15)
print(g)
dev.off()

##roc curve
library(pROC)

rocs <- lapply(columns, function(x) roc(table$perturb, 1/table[[x]], x))
rocs <- lapply(1:length(rocs), function(x) data.frame(sensitivity = rocs[[x]]$sensitivities, specificity = rocs[[x]]$specificities, Method = names[[x]]))
rocs<-do.call(rbind, rocs)

tiff('./roc.tiff', height= 3000, width= 4000, res=500)
g <- ggplot(rocs, aes(1-specificity, sensitivity, colour=Method)) + geom_path() + scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue', 'purple', 'purple4')) + theme_grey(15)
print(g)
dev.off()


