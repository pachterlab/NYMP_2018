
source('~/aggregationDE/R/aggregation.R')
library(dplyr)
library(lmtest)
library(tximport)

base_dir <- '~/single/simulations_genes'
de_dir <- file.path(base_dir, 'de')
folders <- c('rsem_quant_nonperturb', 'rsem_quant_perturb')
files <- sapply(folders, function(folder) sapply(1:105, function(x) file.path(base_dir, folder, paste0(x), 'abundance.h5')))
files <- as.vector(files)
print(length(files))
tximp <- tximport(files, type='kallisto', txOut=TRUE)
saveRDS(tximp, file.path(base_dir, 'tx_tximport.rds'))
counts <- tximp$counts
colnames(counts) <- 1:210

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
colnames(t2g) <- c('transcripts', 'genes')

LR <- function(counts, cutoff = 189)
{
	if(is.vector(counts))
	{
		if(sum(counts  == 0) > cutoff)
		{
			return(NA)
		}
	}
	else{
		zeros <- apply(counts, 2, function(x) sum(x ==0) > cutoff)
		if(sum(zeros) == ncol(counts))
		{
			return(1)
		}
		counts <- counts[,!zeros]
	}
	factors <- as.factor(c(rep(0, 105), rep(1,105)))
	table <- data.frame(counts, factors)
	full <-  glm(factors~., data=table, family=binomial)
	null <-  glm(factors~1, data=table, family=binomial)
	lrtest <- lrtest(full, null)
	lrtest$Pr[2]
}

LR_fit <- function(counts, genes)
{
	table <- data.frame(genes= genes, index = 1:length(genes))
	indices <- table %>% group_by(genes) %>% summarise(indices = list(index))
	pvalues <- sapply(indices$indices, function(x) {LR(counts[,x])})
	n_transcripts <- sapply(indices$indices, length)
	data.frame(pvalues, genes =  indices$genes, n_transcripts)
}

stopifnot(all(t2g$transcripts == rownames(counts)))
counts <- t(counts)
LR_counts <- LR_fit(counts, t2g$genes)
saveRDS(LR_counts, file.path(de_dir, 'LR_tx.rds'))

tpms <- tximp$abundance
colnames(tpms) <- 1:210
stopifnot(all(t2g$transcripts == rownames(tpms)))
LR_tpms <- LR_fit(t(tpms), t2g$genes)
saveRDS(LR_tpms, file.path(de_dir, 'LR_tpms.rds'))

library(DESeq2)
deseq_counts <- tximp$counts
colnames(deseq_counts) <- 1:210
rownames(deseq_counts) <- as.character(rownames(deseq_counts))
deseq_counts <- round(deseq_counts)
condition <- c(rep(0,105), rep(1,105))
condition <- factor(condition)
colData <- data.frame(row.names=colnames(deseq_counts), condition)
dds <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = colData, design = ~condition)
dds <- estimateSizeFactors(dds)
sizeFactors <- dds$sizeFactor
norm_deseq_counts <- sweep(deseq_counts, 2, sizeFactors, '/')
LR_deseq <- LR_fit(t(norm_deseq_counts), t2g$genes)
saveRDS(LR_deseq, file.path(de_dir, 'LR_deseq.rds'))

tests <- apply(counts, 2, function(x) ks.test(x[1:105], x[106:210]))
pvals <- sapply(tests, function(x) x$p.value)
table <- data.frame(pvals = pvals, genes = t2g$genes, weights = colMeans(counts))
results <- table %>% group_by(genes) %>% summarise(pval = lancaster(pvals, weights))
saveRDS(results, file.path(de_dir, 'lan_tx.rds'))


