source('../R/LRtx.R')

tpms <- read.csv('/home/vasilis/DE/Paper_Jun2018/Petropoulos16/EMTAB3929_TPM_matrix.csv', header=TRUE, stringsAsFactors=FALSE)

rownames(tpms) <- tpms[,1]
tpms <- tpms[,-1]

read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', stringsAsFactors=FALSE) -> t2g
names(t2g) <- c('transcripts','genes')


tpms <- readRDS('./embryonic_tpms_matrix.rds')
labels <- read.csv('/home/vasilis/DE/Paper_Jun2018/Petropoulos16/EMTAB3929_labels.csv')
day3 <- which(labels[,2] == 'embryonic day 3')
day4 <- which(labels[,2] == 'embryonic day 4')
day3_4 <- c(day3, day4)
day <- c(rep('3', length(day3)), rep('4', length(day4)))
tpms <- tpms[,day3_4]

coldata <- data.frame(day = day, sample = colnames(tpms))

library('biomaRt')
ensembl = useMart('ensembl', dataset='hsapiens_gene_ensembl')
t2g <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'), mart=ensembl)

transcripts <- rownames(tpms)
transcripts <- sapply(transcripts, function(x) gsub('\\..*', '', x))
rownames(tpms) <- transcripts
tpms <- tpms[rownames(tpms) %in%  t2g$ensembl_transcript_id,]
t2g <- t2g[t2g$ensembl_transcript_id %in% rownames(tpms),]

tpms <- tpms[match(t2g$ensembl_transcript_id,rownames(tpms)), ]
all.equal(t2g$ensembl_transcript_id, rownames(tpms))

genes <- unique(t2g$ensembl_gene_id)
genetpms <- lapply(genes, function(g) {
	i <- which(t2g$ensembl_gene_id == g)
	tpms_g <- tpms[i,]
	if(!is.vector(tpms_g)) {tpms_g <- colSums(tpms_g)}
	tpms_g
})
genetpms <- do.call(rbind, genetpms)
rownames(genetpms) <- genes


library(Seurat)

nbt = CreateSeuratObject(raw.data=round(genetpms), min.cells=3, min.genes=200, project='Embryo')
nbt@ident <- as.factor(coldata[,1])
names(nbt@ident) <- coldata[,2]
deseq2 <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='DESeq2')
saveRDS(deseq2, './deseq2.rds')

mast <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='MAST')
saveRDS(mast, './mast.rds')

tobit <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='tobit')
saveRDS(tobit, './tobit.rds')

source('./LRtx.R')
LR <- LR_fit(tpms, t2g$ensembl_gene_id, ~coldata$day)


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

