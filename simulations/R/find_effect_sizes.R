#Use DESeq2 to find effect sizes for simulation based on effect sizes
library(DESeq2)
library('simulation_framework.R')
counts <- readRDS('../hsmm_counts')
c1 <- lapply(counts, function(x) x[[1]])
c2 <- lapply(counts, function(x) x[[2]])

c1 <- map_sims(c1)
c2 <- map_sims(c2)

c1 <- do.call(cbind, c1)
c2  <- do.call(cbind, c2)

counts <- cbind(c1, c2)
counts <- round(counts)

sampleNames <- c(paste('C1', 1:105), paste('C2', 1:96))
colnames(counts) <- sampleNames

design <- data.frame(row.names = colnames(counts), condition = c(rep('c1', 105), rep('c2', 96)))

DESeqDataSetFromMatrix(countData= counts, colData = design, design=~condition) -> dds

dds <- DESeq(dds)
res <- results(dds)
saveRDS(res, './deseq_hsmm_results.rds')

