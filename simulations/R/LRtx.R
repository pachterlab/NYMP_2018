library(lmtest)

LR <- function(counts, cluster, filter = .9)
{
	cutoff = filter*dim(counts)[[1]]
	if(is.vector(counts))
	{
		if(sum(counts  == 0) > cutoff)
		{
			return(NA)
		}
	}
	else{
		zeros <- apply(counts, 2, function(x) sum(x ==0) > cutoff)
		if(sum(zeros) == ncol(counts))
		{
			return(NA)
		}
		counts <- counts[,!zeros]
	}
	table <- data.frame(counts, cluster)
	full <-  glm(factors~., data=table, family=binomial)
	null <-  glm(factors~1, data=table, family=binomial)
	lrtest <- lrtest(full, null)
	lrtest$Pr[2]
}

#have to generalize
#genes = to what gene does each col of counts belongs to
#cluster = to what cluster each row of counts belongs to
LR_fit <- function(counts, genes, cluster)
{
	table <- data.frame(genes= genes, index = 1:length(genes))
	indices <- table %>% group_by(genes) %>% summarise(indices = list(index))
	pvalues <- sapply(indices$indices, function(x) {LR(counts[,x], cluster)})
	data.frame(pvalues, genes =  indices$genes)
}


