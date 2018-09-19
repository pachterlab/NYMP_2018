

# Supp Fig 5b and 5c
benchmarks <- readRDS('~/single/splat2/benchmark2/benchmark.rds')

library(ggplot2)

levels <- c('DESeq2', 'Monocle - Tobit', 'MAST', 'SCDE', 'log reg - transcripts', 'log reg - gene counts')

tiff('./benchmark.tiff', height = 5000, width = 3500, res = 500)
g <- ggplot(benchmarks, aes(x=factor(Method, levels), y = (Time/60), colour=factor(Method, levels)))+ geom_point(size=1.3) +  scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue', 'purple')) + guides(colour=FALSE)
g <- g+xlab('Method') + ylab('Total CPU Time (min)')
g <- g + theme_bw(base_size=25) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))
g
dev.off()


tiff('./realtime.tiff', height = 5000, width = 3500, res = 500)
g <- ggplot(benchmarks, aes(x=factor(Method, levels), y = (elapsed/60), colour=factor(Method, levels)))+ geom_point(size=1.3) +  scale_colour_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", 'blue', 'purple')) + guides(colour=FALSE)
g <- g + xlab('Method') + ylab('Real Elapsed Time (min)')
g <- g + theme_bw(base_size=25) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))
g
dev.off()


