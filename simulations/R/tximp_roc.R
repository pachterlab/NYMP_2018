source('./roc_curve.R')
library(ggplot2)

base_dir <- '~/single/6.3_sims_exp'
de_dir <- file.path(base_dir, 'de')

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
colnames(t2g) <- c('transcripts', 'genes')

deseq_1 <- readRDS(file.path(de_dir, 'deseq2_1.rds'))
deseq_2 <- readRDS(file.path(de_dir, 'deseq2_2.rds'))
deseq_3<- readRDS(file.path(de_dir, 'deseq2_3.rds'))
LR_tx <- readRDS(file.path(de_dir, 'LR_deseq.rds'))

deseq_1$genes <- rownames(deseq_1)
deseq_2$genes <- rownames(deseq_2)
deseq_3$genes<- rownames(deseq_3)

deseq_1$deseq_1 <- deseq_1$p_val
deseq_2$deseq_2 <- deseq_2$p_val
deseq_3$deseq_3 <- deseq_3$p_val
LR_tx$LR_tx <- LR_tx$pval

merge(deseq_1, deseq_2, by = 'genes', all=TRUE) -> table
merge(deseq_3, table, by ='genes', all=TRUE) -> table
merge(LR_tx, table, by='genes', all=TRUE) -> table

table <- select(table, genes, deseq_1, deseq_2, deseq_3, LR_tx)
perturb <- read.table(file.path(base_dir, 'perturb.tsv'), header=TRUE)

t2g$perturbed <- t2g$transcripts %in% perturb$perturbed_transcripts
perturb <- t2g %>% group_by(genes) %>% summarise(perturb = any(perturbed))
merge(perturb, table, by='genes', all = TRUE) -> table
saveRDS(table, file.path(de_dir, 'tximp_comparison.rds'))

calculate_fdr(table$perturb, table$deseq_1, 'DESeq2 Counts') -> deseq_1
calculate_fdr(table$perturb, table$deseq_2, 'DESeq2 Scaled TPM') -> deseq_2
calculate_fdr(table$perturb, table$deseq_3, 'DESeq2 Length Scaled TPM') -> deseq_3
calculate_fdr(table$perturb, table$LR_tx, 'log reg - transcripts') -> LR_tx

fdrs <- list(deseq_1, deseq_2, deseq_3, LR_tx)
names <- c('DESeq2 - tximport Counts', 'DESeq2 - tximport Scaled TPM', 'DESeq2 - tximport Length Scaled TPM', 'log reg - transcripts')
fdrs <- lapply(1:length(fdrs), function(x) data.frame(sensitivity=fdrs[[x]]$sen, FDR = fdrs[[x]]$fdr, Method = names[[x]])) 
fdrs <- do.call(rbind, fdrs)

p <- ggplot(data = fdrs, aes(x = FDR , y = sensitivity, colour = Method)) + geom_path() + scale_colour_manual(values=c("#F8766D", "#CC0066", "#993366", 'blue'))
tiff(file.path(base_dir, 'tximp_fdrs.png'), height=3000, width=3500, res=600)
p <- p+theme_grey(10)
print(p)
dev.off()

p <- p + coord_cartesian(xlim = c(0, .25), ylim = c(0, .75), expand = TRUE)
png(file.path(base_dir, 'tximp_fdrs2.png'))
print(p)
dev.off()

