
#generate Supp Figi 4 that compares mean variance and dropout rate between simulations and real data
meanVarDF <- readRDS('./mean_variance_comparison.rds')

library(ggplot2)

newDF <- data.frame(means=c(df$exp_means, df$sim_means), var=c(df$exp_var, df$sim_var), data=c(rep('experiment', 173259), rep('simulation', 173259)), zeros=c(df$expzeros, df$simzeros))
newDF$logmeans <- log(newDF$means + 1)
newDF$logvar <- log(newDF$var+1)

g  <- ggplot(newDF, aes(x=logmeans, y=logvar, colour=data)) + geom_point(alpha=.01)
g  <- g + xlab('Log(Mean+1)') + ylab('Log(Var+1)')
png('./meanVar.png')
print(g)
dev.off()

g <- ggplot(newDF, aes(x=logmeans, y=zeros/105, colour=data)) + geom_point(alpha=.01)
g  <- g + xlab('Log(Mean+1)') + ylab('Percentage Zeros')
png('./dropout.png')
print(g)
dev.off()


diffDF <- data.frame(diffMeans = df$lsimmeans - df$lexpmean, diffVar=var(df$lsimmeans) - var(df$lexpmean), diffZeros = (df$simzeros - df$expzeros)/105)

library(reshape2)
diffDF  <- melt(diffDF)

g <- ggplot(diffZeros, aes(x=variable, y=value)) + geom_boxplot()
png('./diff.png')
print(g)
dev.off()

diffDF <- data.frame(diffMeans = df$lsimmeans - df$lexpmean, diffVar=var(df$lsimmeans) - var(df$lexpmean))


