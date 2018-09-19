
library(tximport)
library(Seurat)

args <- commandArgs(TRUE)
base_dir = args[[1]] 
#base_dir <- '/home/lynnyi/single/6.3_sims_exp/'
de_dir = file.path(base_dir, 'de')

if(file.exists(file.path(base_dir, 'tximport.rds'))){
	tximp <-readRDS(file.path(base_dir, 'tximport.rds'))
} else{
	folders <- c('rsem_quant_nonperturb', 'rsem_quant_perturb')
	folders <- sapply(folders, function(folder) file.path(base_dir, folder))
	files <- sapply(folders, function(folder) sapply(1:105, function(x) file.path(folder, paste0(x), 'abundance.h5')))
	files <- as.vector(files)
	print(files)
	print(length(files))
	t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts')
	tximp <- tximport(files, type='kallisto', tx2gene=t2g, txOut=FALSE)
	saveRDS(tximp, file.path(base_dir, 'tx_tximport.rds'))
}

tpms <- tximp$abundance
colnames(tpms) <- 1:210
counts <- round(tximp$counts)
colnames(counts) <- 1:210

nbt = CreateSeuratObject(raw.data = counts, min.cells=0, min.genes=0, project='Sim')
nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
names(nbt@ident) <- rep(1:210)
deseq2 <- FindMarkers(nbt, ident.1= 0, min.pct = .1, test.use = 'DESeq2')
saveRDS(deseq2, file.path(de_dir, 'deseq2.rds'))

#nbt = CreateSeuratObject(raw.data = tpms, min.cells=0, min.genes=0, project='Sim')
#nbt <- NormalizeData(nbt)
#nbt <- ScaleData(nbt)
#nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
#names(nbt@ident) <- rep(1:210)
#nbt.markers <- FindAllMarkers(nbt, min.pct = .1)
#saveRDS(nbt.markers, file.path(de_dir, 'seurat.rds'))
#
#nbt = CreateSeuratObject(raw.data = tpms, min.cells=0, min.genes=0, project='Sim')
#nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
#names(nbt@ident) <- rep(1:210)
#time_mast <- system.time(mast <- FindMarkers(nbt, ident.1 = 0, min.pct=.1,  test.use = 'MAST'))
#saveRDS(mast, file.path(de_dir, 'mast.rds'))
#

tobit <- FindMarkers(nbt, ident.1 = 0, min.pct=.1,  test.use = 'tobit')
saveRDS(tobit, file.path(de_dir, 'tobit.rds'))

nbt = CreateSeuratObject(raw.data = log(tpms+1), min.cells=0, min.genes=0, project='Sim')
nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
names(nbt@ident) <- rep(1:210)
time_mast <- system.time(mast <- FindMarkers(nbt, ident.1 = 0, min.pct=.1,  test.use = 'MAST'))
saveRDS(mast, file.path(de_dir, 'mast_log.rds'))


