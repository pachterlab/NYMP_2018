source('./roc_curve.R')
library(ggplot2)

base_dir <- '~/single/simulations_genes'
de_dir <- file.path(base_dir, 'de')
scde_dir <- '/home/scde/simulations_genes/'

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
colnames(t2g) <- c('transcripts', 'genes')

deseq <- readRDS(file.path(de_dir, 'deseq2.rds'))
mast <- readRDS(file.path(de_dir, 'mast.rds'))
lan <- readRDS(file.path(de_dir, 'lan.rds'))
LR <- readRDS(file.path(de_dir, 'LR.rds'))
LR_tx <- readRDS(file.path(de_dir, 'LR_tx.rds'))
scde <- readRDS(file.path(scde_dir, 'scde_DE.rds'))
tobit <- readRDS(file.path(de_dir, 'tobit.rds'))
LR_deseq <- readRDS(file.path(de_dir, 'LR_deseq.rds'))
LR_tpms <- readRDS(file.path(de_dir, 'LR_tpms.rds'))


deseq$genes <- rownames(deseq)
mast$genes <- rownames(mast)
scde$genes <- rownames(scde)
tobit$genes<- rownames(tobit)

deseq$deseq <- deseq$p_val
mast$mast <- mast$p_val
tobit$tobit <- tobit$p_val
lan$lan <- lan$pval
LR$LR <- LR$pval
LR_tx$LR_tx <- LR_tx$pval
scde$scde <- scde$pval
LR_deseq$LR_deseq <- LR_deseq$pval
LR_tpms$LR_tpms <- LR_tpms$pval
LR_deseq <- select(LR_deseq, genes, LR_deseq)
LR_tpms <- select(LR_tpms, genes, LR_tpms)

merge(deseq, mast, by = 'genes', all=TRUE) -> table
merge(table, scde, by='genes', all=TRUE) -> table
merge(table, LR, by='genes', all=TRUE) -> table
merge(lan, table, by='genes', all=TRUE) -> table
merge(LR_tx, table, by='genes', all=TRUE) -> table
merge(tobit, table, by='genes', all=TRUE) -> table
merge(table, LR_deseq, by='genes', all=TRUE) -> table
merge(table, LR_tpms, by='genes', all=TRUE) -> table


table <- select(table, genes, mast, deseq, scde, lan, LR, LR_tx, n_tccs, n_transcripts, tobit, LR_deseq, LR_tpms)
#table <- select(table, genes, mast, deseq, scde, lan, LR, LR_tx, n_tccs, n_transcripts, LR_deseq, LR_tpms)

table$lr70 <- table$LR
table$lr70[table$n_tccs > 70] <- table$lan[table$n_tccs > 70]

#table$mast <- p.adjust(table$mast, 'BH')
#table$deseq <- p.adjust(table$deseq, 'BH')
#table$scde <- p.adjust(table$scde, 'BH')
#table$tobit <- p.adjust(table$tobit, 'BH')
#table$LR_tx <- p.adjust(table$LR_tx, 'BH')
#table$lr70 <- p.adjust(table$lr70, 'BH')

read.table(file.path(base_dir, 'perturb.tsv'), header=TRUE) -> perturb
t2g$perturbed <- t2g$transcripts %in% perturb$perturbed_transcripts
perturb <- t2g %>% group_by(genes) %>% summarise(perturb = any(perturbed))
merge(perturb, table, by='genes', all = TRUE) -> table
saveRDS(table, file.path(de_dir, 'all_results.rds'))

calculate_fdr(table$perturb, table$mast, 'mast') -> mast
calculate_fdr(table$perturb, table$deseq, 'deseq') -> deseq
calculate_fdr(table$perturb, table$scde, 'scde') -> scde
calculate_fdr(table$perturb, table$tobit, 'tobit') -> tobit

LR_tccs <- calculate_fdr(table$perturb, table$lr70, 'LR_tccs')
LR_tx <- calculate_fdr(table$perturb, table$LR_tx, 'LR_tx')
LR_deseq <- calculate_fdr(table$perturb, table$LR_deseq, 'LR_deseq')
LR_tpms <- calculate_fdr(table$perturb,table$LR_tpms, 'LR_tpms')

plot_fdrs(mast,tobit,  deseq, scde, LR_tccs, LR_tx, file='~/fdrs.png', title = 'test')

fdrs <- list(deseq, tobit, mast, scde, LR_tx)
names <- c('DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE', 'LogR')
fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]])) 
fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue'))
png(file.path(base_dir, 'fdrs.png'))
print(p)
dev.off()

p <- p + coord_cartesian(xlim = c(0, .25), ylim = c(0, .75), expand = TRUE)
png(file.path(base_dir, 'fdrs2.png'))
print(p)
dev.off()

fdrs <- list(LR_tx, LR_tccs, deseq, tobit, mast, scde)
names <- c('LogR - Transcripts', 'LogR - TCCs', 'DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE') 

fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]]))
fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values = c('blue', 'red', 'grey', 'grey', 'grey', 'grey'))

png(file.path(base_dir, 'fdrs3.png'))
print(p)
dev.off()


fdrs <- list(LR_tx, LR_deseq, LR_tpms, deseq, tobit, mast, scde)
names <- c('LogR - Counts', 'LogR - DESeq2 Norm', 'LogR - TPM Norm','DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE')

fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]]))

fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values = c('blue', 'darkcyan', 'blueviolet', 'grey', 'grey', 'grey', 'grey'))

png(file.path(base_dir, 'fdrs4.png'))
print(p)
dev.off()

p <- p + coord_cartesian(xlim = c(0, .25), ylim = c(0, .75), expand = TRUE)
png(file.path(base_dir, 'fdrs5.png'))
print(p)
dev.off()
