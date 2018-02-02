
library(dplyr)

calculate_fdr <- function(true_DE, pvalues, title)
{
	df <- data.frame(true_DE = true_DE, pvalues=pvalues)
	df <- df[order(df$pvalues), ]
	total_positive <- sum(true_DE)
	df <- df %>% group_by(pvalues) %>% summarise(p = sum(true_DE), n = sum(!true_DE), l=n())
	print(head(df))
	fdr <- sapply(seq(df$pvalues), function(i) sum(df$n[1:i]) / sum(df$l[1:i]))
	sensitivity <- sapply(seq(df$pvalues), function(i) sum(df$p[1:i])/total_positive)
	n_found <- sapply(seq(df$pvalues), function(i) sum(df$p[1:i]))
	n_de <- sapply(seq(df$pvalues), function(i) sum(df$l[1:i]))
	five <- min(which(df$pvalues > .05)) -1 
	ten <- min(which(df$pvalues > .1)) - 1
	twenty <- min(which(df$pvalues > .2)) -1
	
	list(fdr = fdr, sensitivity = sensitivity,
		n_found = n_found,
		n_de = n_de,
		pvalues = df$pvalues,
		five=five, ten=ten, twenty=twenty, title=title)
}


