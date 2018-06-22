library(tximport)
library(Seurat)

base_dir = '~/single/6.3_sims_exp/'
de_dir = file.path(base_dir, 'de')

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts')
folders <- c('rsem_quant_nonperturb', 'rsem_quant_perturb')
folders <- sapply(folders, function(folder) file.path(base_dir, folder))

files <- sapply(folders, function(folder) sapply(1:105, function(x) file.path(folder, paste0(x), 'abundance.h5')))

files <- as.vector(files)
print(files)
print(length(files))
tximp1 <- tximport(files, type = 'kallisto', tx2gene = t2g, countsFromAbundance='no')
tximp2 <- tximport(files, type='kallisto', tx2gene=t2g, countsFromAbundance='scaledTPM')
tximp3 <- tximport(files, type='kallisto', tx2gene=t2g, countsFromAbundance='lengthScaledTPM')
tximps <- list(tximp1, tximp2, tximp3)
sapply(1:3, function(x) {
	counts <- round(tximps[[x]]$counts)
	colnames(counts) <- 1:210
	nbt = CreateSeuratObject(raw.data = counts, min.cells=3, min.genes=200, project='Sim')
	nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
	names(nbt@ident) <- rep(1:210)
	deseq2 <- FindMarkers(nbt, ident.1= 0, min.pct = .25, test.use = 'DESeq2')
	saveRDS(deseq2, file.path(de_dir, paste0('deseq2_', x, '.rds')))
})

