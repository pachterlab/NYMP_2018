library(aggregation)

make_tcc_to_gene_map <- function(transcripts, map)
{
        ntccs <- nrow(map)
        tcc2gene <- sapply(1:ntccs, function(i) map_to_gene(i, transcripts, map))
        tcc2gene
}

map_to_gene <- function(i, transcripts, map)
{
        #map tcc i to transcript id
        transcript_indices <- unlist(strsplit(toString(map[[i,2]]), ','))
        transcript_indices <- sapply(transcript_indices, strtoi)
        #have to map from 0 indexed transcript indices to the t2g map
        transcript_indices <- transcript_indices +1
        g <- sapply(transcript_indices, function(x) transcripts[x,2])
        unique_g <- unique(g)

        if(length(unique_g) == 1)
        {
                return(unique_g)
        }
        else
        {
                return('')
        }
}

library(dplyr)
library(lmtest)
library(Matrix)
library(DESeq2)

args <- commandArgs(TRUE)
base_dir <- args[[1]]
#base_dir <- '~/single/6.3_sims_exp/'
de_dir <- file.path(base_dir, 'de')
do_deseq_norm <- TRUE

print('reading matrix')
m <- read.table(file.path(base_dir, 'pseudo_batch', 'matrix.tsv'))
m <- sparseMatrix(i = m$V2 + 1, j = m$V1 +1, x = m$V3)
m <- as.matrix(m)
print(dim(m))

#normalize with DESeq2 for good measure
if(do_deseq_norm)
{
	print('normalizing')
	source('~/single/R/deseq_norm.R')
	m <- t(deseq_norm(t(m)))
}

if(file.exists(file.path(base_dir, 'tcc2gene.rds'))){
	print('reading map')
	map <- readRDS(file.path(base_dir, 'tcc2gene.rds'))
	if(typeof(map)=='list')
		{map <- unlist(map)}
} else{
	ec_map <- read.table(file.path(base_dir, 'pseudo_batch', 'matrix.ec'), stringsAsFactors=FALSE, sep='\t')
	t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
	map <- make_tcc_to_gene_map(t2g, ec_map)
	saveRDS(map, file.path(base_dir, 'tcc2gene.rds'))
}

print('performing KS testing')
tests <- apply(m, 2, function(x) ks.test(x[1:105], x[106:210]))
pvals <- sapply(tests, function(x) x$p.value)
table <- data.frame(pvals = pvals, genes = map, weights = colMeans(m))
results <- table %>% group_by(genes) %>% summarise(pval = lancaster(pvals, weights))
results <- results[-1,] 
saveRDS(results, file.path(de_dir, 'lan_tccs.rds'))

LR <- function(tcc_counts)
{
	factors <- as.factor(c(rep(0, 105), rep(1,105)))
	table <- data.frame(tcc_counts)
	i <- apply(table, 2, function(x) sum(x==0) < 189)
	if(sum(i) == 0)
	{
		return(list(p = NA, n_tccs= 0))
	}
	if(sum(i)==1)
	{
		table <- as.data.frame(table[,i])
	}
	else
	{
		table <- table[,i]
	}
	table$factors <- factors
	full <-  glm(factors~., data=table, family=binomial)
	null <-  glm(factors~1, data=table, family=binomial)
	lrtest <- lrtest(full, null)
	return(list(p = lrtest$Pr[2], n_tccs = sum(i)))
}

LR_fit <- function(counts, tcc_map)
{
	table <- data.frame(genes= tcc_map, index = 1:length(tcc_map))
	indices <- table %>% group_by(genes) %>% summarise(indices = list(index))
	indices <- indices[-1,]
	results <- lapply(indices$indices, function(x) {LR(counts[,x])})
	pvalues <- sapply(results, function(x) x$p)
	n_tccs <- sapply(results, function(x) x$n_tccs)
	data.frame(pvalues, genes =  indices$genes, n_tccs)
}

print('performing log r')
LR <- LR_fit(m, map)
saveRDS(LR, file.path(de_dir, 'LR_tccs.rds'))

