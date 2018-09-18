
#Load and parse data ----------------------------------------------------
library(MultiAssayExperiment)
data <- readRDS('/home/lynnyi/NYMP_2018/embryo/EMTAB3929.rds')
gene <- experiments(data)[['gene']]
tpms <- assays(gene)[['TPM']][,1:271]
counts <- assays(gene)[['count']][,1:271]
day <- c(rep(3, 81), rep(4, 190))

# Use Seurat to call DESeq2 ---------------------------------------------
library(Seurat)
nbt = CreateSeuratObject(raw.data=round(counts), min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(counts) 
deseq2 <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='DESeq2')
saveRDS(deseq2, './deseq2.rds')

# Use Seurat to call MAST -----------------------------------------------
nbt = CreateSeuratObject(raw.data=log(tpms+1), min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(tpms)
mast <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='MAST')
saveRDS(mast, './mast.rds')

# Use Seurat to run Tobit model --------------------------------------------
nbt = CreateSeuratObject(raw.data=tpms, min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(tpms)
tobit <- FindMarkers(nbt, ident.1=3, min.pct=.1, test.use='tobit')
saveRDS(tobit, './tobit.rds')

