source('../R/LRtx.R')

# Load and parse data ------------------------------------------------
data <- readRDS('/home/lynnyi/NYMP_2018/embryo/EMTAB3929.rds')
tx <- experiments(data)[['tx']]
tx_counts <- assays(tx)[['count']][,1:271]

# Get transcript to gene table --------------------------------------
library('biomaRt')
ensembl = useMart('ensembl', dataset='hsapiens_gene_ensembl')
t2g <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
	mart=ensembl)

rownames(tx_counts) <- gsub('\\..*', '', rownames(tx_counts))
tx_counts <- tx_counts[rownames(tx_counts) %in%  t2g$ensembl_transcript_id,]
t2g <- t2g[t2g$ensembl_transcript_id %in% rownames(tx_counts),]

tx_counts <- tx_counts[match(t2g$ensembl_transcript_id,rownames(tx_counts)), ]
all.equal(t2g$ensembl_transcript_id, rownames(tx_counts))

# Perform logistic regression on transcript quantifications -----------
source('~/NYMP_2018/simulations/R/LRtx.R')
source('~/NYMP_2018/simulations/R/deseq_norm.R')
tx_counts_norm <- deseq_norm(tx_counts)
day <- as.factor(c(rep(0, 81), rep(1, 190)))
LR <- LR_fit(t(tx_counts_norm), t2g$ensembl_gene_id, day, filter=.9)


# Load other methods' results ----------------------------------------------
scde <- readRDS('./scde.rds')

results <- list(scde, mast, tobit, deseq2)
names <- list('scde', 'mast', 'tobit', 'deseq2')
results <- lapply(1:length(results), function(x){
	results[[x]]$genes <- rownames(results[[x]])
	results[[x]][[names[[x]]]] <- results[[x]]$p_val
	if(!names[[x]] %in% names(results[[x]]))
		results[[x]][[names[[x]]]] <- results[[x]]$pval
	results[[x]]})


LR$LR <- LR$pvalues
table <- merge(LR, results[[1]], by='genes', all=TRUE)
for(i in 2:length(results))
{
	table <- merge(table, results[[i]], by='genes', all=TRUE)
}
rownames(table) <- table$genes
pvals <- table[,which(colnames(table) %in% c(names, 'LR'))]
names(pvals) <- c('logistic regression', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
adjusted <- apply(pvals, 2, function(x) p.adjust(x))
sig_genes <- apply(pvals, 2, function(x) rownames(pvals)[which(x < .05)])
sig_genes_adjusted <- apply(adjusted, 2, function(x) rownames(pvals)[which(x<.05)])

bon <- adjusted
bon[,1] <- p.adjust(pvals[,1], method= 'bonferroni')
sig_genes_bon <- apply(bon, 2, function(x) rownames(pvals)[which(x<.05)])

library(UpSetR)
tiff('./upset.png', height=3000, width=3000, res = 500)
upset(fromList(sig_genes_adjusted), order.by='freq')
dev.off()

tiff('./upset.png', height=3000, width=3000, res = 500)
upset(fromList(sig_genes_bon), order.by='freq')
dev.off()

# lots of deprecated genes that do not exist in existing annotation
scde <- scde[rownames(scde) %in% LR$genes,]
tobit <- tobit[rownames(tobit) %in% LR$genes,]
deseq <- deseq[rownames(deseq) %in% LR$genes,]
mast <- mast[rownames(mast) %in% LR$genes,]

LR$p_val_adj <- p.adjust(LR$pvalues, 'BH')
scde$p_val_adj <- p.adjust(scde$pvalue, 'BH')
l <- list(LR, scde, tobit, deseq, mast)
test <- lapply(l, function(x) rownames(x[order(x$p_val_adj),])[1:3000])
names(test) <- c('log reg', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
tiff('./upset.tiff', height=3000, width = 4000, res= 500)
upset(fromList(test), order.by='freq', text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()

unique <- setdiff(test[[1]], test[[2]])
unique <- setdiff(unique, test[[3]])
unique <- setdiff(unique, test[[4]])
unique <- setdiff(unique, test[[5]])

plot_gene <- function(gene_name)
{
	i <- which(t2g$ensembl_gene_id == gene_name)
	tx <- t2g$ensembl_transcript_id[i]
	print(length(tx))
	tx_counts <- tx_counts[rownames(tx_counts) %in% tx,]
	tx_counts <- log(tx_counts+1)
	tx_counts <- melt(tx_counts)
	names(tx_counts) <- c('transcript', 'day', 'counts')
	tx_counts$day <- as.factor(tx_counts$day)
	g <- ggplot(tx_counts, aes(transcript, counts, fill=day)) + geom_boxplot(alpha=.2) + ylab('log(counts+1)')
	g <- g+theme(axis.text.x = element_text(angle = 90, hjust = 1))
	g <- g+ggtitle(gene_name)
	tiff(paste0(gene_name,'.tiff'))
	print(g)
	dev.off()
}


