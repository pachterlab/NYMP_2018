
library(MultiAssayExperiment)

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
source('../LRtx.R')
source('~/NYMP_2018/simulations/RSEM/R/deseq_norm.R')
tx_counts_norm <- deseq_norm(tx_counts)
day <- as.factor(c(rep(0, 81), rep(1, 190)))
LR <- LR_fit(t(tx_counts_norm), t2g$ensembl_gene_id, day, .9)
saveRDS('./LR.rds')

