#Perform logistic regression on gene counts for simulations

library(tximport)
library(lmtest)
library(dplyr)
library(DESeq2)

args <- commandArgs(TRUE)
base_dir <- args[[1]]
#base_dir <- '/home/lynnyi/single/5.22_sims_genes/'
de_dir <- file.path(base_dir, 'de')

tximp <- readRDS(file.path(base_dir, 'tximport.rds'))
counts <- tximp$counts

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
                        return(NA)
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

counts <- round(counts)
condition <- factor(c(rep(0, 105), rep(1, 105)))
colData <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors <- dds$sizeFactor
counts <- sweep(counts, 2, sizeFactors, '/')

LR_uni_time <- system.time(pvals <- apply(counts, 1, LR))
saveRDS(LR_uni_time, file.path(de_dir, 'LR_uni_time2.rds'))
results <- data.frame(genes = rownames(counts), pvals = pvals)
saveRDS(results, file.path(de_dir, 'LR_univariate_deseq2.rds'))

