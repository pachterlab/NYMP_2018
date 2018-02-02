
library(aggregation)
library(dplyr)
library(lmtest)

base_dir = '~/single/sims_genes/'
de_dir = file.path(base_dir, 'de')
m <- scan(file.path(base_dir, 'pseudo_batch', 'matrix.tsv'), skip = 1, what=numeric(), sep='\t')
m <- matrix(m, nrow=211)
m <- m[-1,]
print(dim(m))

ec_map <- read.table(file.path(base_dir, 'pseudo_batch', 'matrix.ec'), stringsAsFactors=FALSE, sep='\t')
t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE, sep='\t')
print(typeof(t2g))
map <- make_tcc_to_gene_map(t2g, ec_map)

tests <- apply(m, 2, function(x) ks.test(x[1:105], x[106:210]))
pvals <- sapply(tests, function(x) x$p.value)
table <- data.frame(pvals = pvals, genes = map, weights = colMeans(m))
results <- table %>% group_by(genes) %>% summarise(pval = lancaster(pvals, weights))
results <- results[-1,] 
saveRDS(results, file.path(de_dir, 'lan.rds'))

LR <- function(tcc_counts)
{
	factors <- as.factor(c(rep(0, 105), rep(1,105)))
	table <- data.frame(tcc_counts)
	table$factors <- factors
	full <-  glm(factors~., data=table, family=binomial)
	null <-  glm(factors~1, data=table, family=binomial)
	lrtest <- lrtest(full, null)
	lrtest$Pr[2]
}

LR_fit <- function(counts, tcc_map)
{
	table <- data.frame(genes= tcc_map, index = 1:length(tcc_map))
	indices <- table %>% group_by(genes) %>% summarise(indices = list(index))
	indices <- indices[-1,]
	pvalues <- sapply(indices$indices, function(x) {LR(counts[,x])})
	n_tccs <- sapply(indices$indices, length)
	data.frame(pvalues, genes =  indices$genes, n_tccs)
}


LR <- LR_fit(m, map)
saveRDS(LR, file.path(de_dir, 'LR.rds'))

