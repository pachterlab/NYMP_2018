
library(tximport)
library(Seurat)

base_dir = '~/single/sims_genes/'
de_dir = file.path(base_dir, 'de')

t2g <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts')
folders <- c('rsem_quant_nonperturb', 'rsem_quant_perturb')
folders <- sapply(folders, function(folder) file.path(base_dir, folder))

files <- sapply(folders, function(folder) sapply(1:105, function(x) file.path(folder, paste0(x), 'abundance.h5')))

files <- as.vector(files)
print(files)
print(length(files))
tximp <- tximport(files, type='kallisto', tx2gene=t2g)
saveRDS(tximp, file.path(base_dir, 'tximport.rds'))
tpms <- tximp$abundance
colnames(tpms) <- 1:210
counts <- round(tximp$counts)
colnames(counts) <- 1:210

nbt = CreateSeuratObject(raw.data = counts, min.cells=3, min.genes=200, project='Sim')
nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
names(nbt@ident) <- rep(1:210)
deseq2 <- FindMarkers(nbt, ident.1= 0, min.pct = .25, test.use = 'DESeq2')
saveRDS(deseq2, file.path(de_dir, 'deseq2.rds'))

#nbt <- RunPCA(object = nbt, pc.genes = rownames(tpms), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#nbt <- FindClusters(object = nbt , reduction.type = "pca", dims.use = 1:10,  resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
#PrintFindClustersParams(nbt)


nbt = CreateSeuratObject(raw.data = tpms, min.cells=3, min.genes=200, project='Sim')
nbt <- NormalizeData(nbt)
nbt <- ScaleData(nbt)
nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
names(nbt@ident) <- rep(1:210)
nbt.markers <- FindAllMarkers(nbt, min.pct = .25, thresh.use=0.25)
saveRDS(nbt.markers, file.path(de_dir, 'seurat.rds'))

nbt = CreateSeuratObject(raw.data = tpms, min.cells=3, min.genes=200, project='Sim')
nbt@ident <- as.factor(c(rep(0,105), rep(1,105)))
names(nbt@ident) <- rep(1:210)
mast <- FindMarkers(nbt, ident.1 = 0, min.pct=.25,  test.use = 'MAST')
saveRDS(mast, file.path(de_dir, 'mast.rds'))

tobit <- FindMarkers(nbt, ident.1 = 0, min.pct=.25,  test.use = 'tobit')
saveRDS(tobit, file.path(de_dir, 'tobit.rds'))

