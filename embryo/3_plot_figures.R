
# Load methods' results ----------------------------------------------
scde <- readRDS('./scde.rds')
mast <- readRDS('./mast.rds')
tobit <- readRDS('./tobit.rds')
deseq2 <- readRDS('./deseq.rds')
LR <- readRDS('./LR.rds')

# Make sure the genes match, since some pulled from conquer and some genes (LR) pulled from ensembl  ----------------------------------------
scde <- scde[rownames(scde) %in% LR$genes,]
tobit <- tobit[rownames(tobit) %in% LR$genes,]
deseq <- deseq[rownames(deseq) %in% LR$genes,]
mast <- mast[rownames(mast) %in% LR$genes,]

# Perform multiple hypothesis adjustment ----------------------------
LR$p_val_adj <- p.adjust(LR$pvalues, 'BH')
scde$p_val_adj <- p.adjust(scde$pvalue, 'BH')

# Create upset plot Supp Fig 6 --------------------------------
# Take top 3000 genes from all methods
l <- list(LR, scde, tobit, deseq, mast)
test <- lapply(l, function(x) rownames(x[order(x$p_val_adj),])[1:3000])
names(test) <- c('log reg', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
tiff('./upset.tiff', height=3000, width = 4000, res= 500)
upset(fromList(test), order.by='freq', text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()

# Find what is unique to logistic regression
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

# Supp Fig 6b
plot_gene('ENSG00000204481')

# Supp Fig 6c
plot_gene('ENSG00000128253')

# Supp Fig 6d
plot_gene('ENSG00000122565')

# Supp Fig 6e
plot_gene('ENSG00000106682')
