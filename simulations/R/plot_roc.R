source('/home/lynnyi/single/R/roc_curve.R')
library(ggplot2)

args <- commandArgs(TRUE)
base_dir <- args[[1]]
#base_dir <- '~/single/6.3_sims_exp/'
de_dir <- file.path(base_dir, 'de')
adjust <- FALSE

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
colnames(t2g) <- c('transcripts', 'genes')

deseq <- readRDS(file.path(de_dir, 'deseq2.rds'))
mast <- readRDS(file.path(de_dir, 'mast_log.rds'))
scde <- readRDS(file.path(de_dir, 'scde.rds'))
tobit <- readRDS(file.path(de_dir, 'tobit.rds'))

#TCCs
lan <- readRDS(file.path(de_dir, 'lan_tccs.rds'))
LR <- readRDS(file.path(de_dir, 'LR_tccs.rds'))
#different normalizations on transcripts
LR_tx <- readRDS(file.path(de_dir, 'LR_tx.rds'))
LR_deseq <- readRDS(file.path(de_dir, 'LR_deseq.rds'))
LR_tpms <- readRDS(file.path(de_dir, 'LR_tpms.rds'))
#LR on genes, DESeq2 normalization
LR_univariate <- readRDS(file.path(de_dir, 'LR_univariate_deseq.rds'))

deseq$genes <- rownames(deseq)
mast$genes <- rownames(mast)
scde$genes <- rownames(scde)
tobit$genes<- rownames(tobit)

#rename for easier merging
deseq$deseq <- deseq$p_val
mast$mast <- mast$p_val
tobit$tobit <- tobit$p_val
lan$lan <- lan$pval
LR$LR <- LR$pval
LR_tx$LR_tx <- LR_tx$pval
scde$scde <- scde$pval
LR_deseq$LR_deseq <- LR_deseq$pval
LR_tpms$LR_tpms <- LR_tpms$pval
LR_univariate$LR_univariate <- LR_univariate$pval



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
merge(table, LR_univariate, by='genes', all=TRUE) -> table

table <- select(table, genes, mast, deseq, scde, lan, LR, LR_tx, n_tccs, n_transcripts, tobit, LR_deseq, LR_tpms, LR_univariate)
n <- 52
table$LR_tccs <-table$LR
table$LR_tccs[table$n_tccs > n] <- NA

if(adjust)
{
	table$mast <- p.adjust(table$mast, 'BH')
	table$deseq <- p.adjust(table$deseq, 'BH')
	table$scde <- p.adjust(table$scde, 'BH')
	table$tobit <- p.adjust(table$tobit, 'BH')
	table$LR_tx <- p.adjust(table$LR_tx, 'BH')
	table$LR_tccs <- p.adjust(table$LR_tccs, 'BH')
	table$LR_univariate <- p.adjust(table$LR_univariate, 'BH')
}

read.table(file.path(base_dir, 'perturb.tsv'), header=TRUE) -> perturb
t2g$perturbed <- t2g$transcripts %in% perturb$perturbed_transcripts
perturb <- t2g %>% group_by(genes) %>% summarise(perturb = any(perturbed))
merge(perturb, table, by='genes', all = TRUE) -> table
saveRDS(table, file.path(de_dir, 'all_results.rds'))

calculate_fdr(table$perturb, table$mast, 'mast') -> mast
calculate_fdr(table$perturb, table$deseq, 'deseq') -> deseq
calculate_fdr(table$perturb, table$scde, 'scde') -> scde
calculate_fdr(table$perturb, table$tobit, 'tobit') -> tobit

LR_tccs <- calculate_fdr(table$perturb, table$LR_tccs, 'LR_tccs')
LR_tx <- calculate_fdr(table$perturb, table$LR_tx, 'LR_tx')
LR_deseq <- calculate_fdr(table$perturb, table$LR_deseq, 'LR_deseq')
LR_tpms <- calculate_fdr(table$perturb,table$LR_tpms, 'LR_tpms')
LR_univariate <- calculate_fdr(table$perturb, table$LR_univariate, 'LR_univariate')

#PLOT ONE
fdrs <- list(deseq, tobit, mast, scde, LR_deseq, LR_univariate, LR_tccs)
names <- c('DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE', 'log reg - transcripts', 'log reg - gene counts', 'log reg - TCCs')
fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]])) 
fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue', 'purple', 'purple4'))
p <- p+theme_gray(15)
tiff(file.path(base_dir, 'fdrs.tiff'), height=3000, width=3500, res=500)
print(p)
dev.off()

#PLOT TWO ZOOM
p <- p + coord_cartesian(xlim = c(0, .25), ylim = c(0, .75), expand = TRUE)
tiff(file.path(base_dir, 'fdrs2.tiff'), height=3000, width=3500, res=500)
p <- p + theme_gray(15)
print(p)
dev.off()


#PLOT FOUR Compare normalizations
fdrs <- list(deseq, tobit, mast, scde, LR_tx, LR_tpms, LR_deseq)
names <- c('DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE', 'log reg - raw counts', 'log reg - TPM Norm', 'log reg - DESeq2 Norm')

fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]]))

fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values = c('grey', 'grey', 'grey', 'grey', 'darkcyan', 'blueviolet', 'blue'))
#try to reverse plot order so show up on top
p <- p + theme_gray(12)

tiff(file.path(base_dir, 'fdrs4.tiff'), height=3000, width=3500, res=600)
print(p)
dev.off()


