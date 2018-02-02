
library(scde)
library(methods)
library(tximport)

source('read.R')

args = commandArgs(TRUE)

nonperturb <- args[[1]]
perturb <- args[[2]]
error_models_path <- args[[3]]
de_path <- args[[4]]

perturb_table_path <- args[[5]]
df_path <- args[[6]]
qqplot_path <- args[[7]]

print('reading kallisto for scde')
dirs <- c(nonperturb, perturb)
print(dirs)
files <- sapply(dirs, function(x) file.path(x, list.files(x), 'abundance.tsv'))
files <- as.vector(unlist(files))
print(length(files))

tx2gene <- read.table('~/transcriptomes/Homo_sapiens.GRCh38.rel79.transcripts', sep='\t')
tximp <- tximport(files, type='kallisto', tx2gene = tx2gene) 

sg <- factor(c(rep(1, 105), rep(2, 105)))
counts <- tximp$counts
counts <- apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
colnames(counts) <- 1:ncol(counts)

print('dim(counts)')
print(dim(counts))
print(head(counts))

print('calculating error models')
error_models <- scde.error.models(counts, sg, n.cores=1, verbose=1, min.count.threshold = 4, min.nonfailed=5, save.model.plots=TRUE)
prior <- scde.expression.prior(models=error_models, counts)

print('writing error models table')
write.table(error_models, error_models_path, sep=',', quote=FALSE)

print('calculating posteriors')
scde_DE <- scde.expression.difference(error_models, counts, prior, sg, verbose=1, n.cores=1, n.randomizations=100)

scde_DE$pvalue <- 2*pnorm(-abs(scde_DE$Z))

print('writing de table')
write.table(scde_DE, de_path, sep=',', quote=FALSE)


source('roc_curve.R')
perturb_table <- read.table(perturb_table_path, header=TRUE, sep='\t')
perturb_transcripts <- perturb_table$perturbed_transcripts
df <- make_data_frame_scde(perturb_transcripts, scde_DE)
save_df(df, df_path)

png(qqplot_path)
qqplot(df$pvalue, 'SCDE QQPlot')
dev.off()


