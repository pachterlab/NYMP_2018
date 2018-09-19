#method of deseq normalization
library(DESeq2)

deseq_norm <- function(m)
{
	ncols <- dim(m)[2]
	condition <- c(0, rep(1, ncols-1))
	m <- round(m)
	colData <- data.frame(row.names=colnames(m), condition = factor(condition))
	dds <- DESeqDataSetFromMatrix(countData=m, colData=colData, design = ~condition)
	dds <- estimateSizeFactors(dds)
	sizeFactors <- dds$sizeFactor
	m <- sweep(m, 2, sizeFactors, '/')
	m	
}
