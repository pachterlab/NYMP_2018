
# need to read without normalizing differently, or maybe by normalizing using DEseq

library(testthat)
library(truncnorm)
library(MASS)
library(lmtest)

source('read.R')

cutoff = 150

fit <- function(x,  model = 'negative binomial', rsem_table)
{
	results <- list()
	if(names(x) != rsem_table$transcript_id)
	{
		print('warning! rsem table order does not match.')
	}
	for(i in 1:length(x))
	{
		result <- tryCatch(
		{
			if (model == 'negative binomial') {
				y <- round(x[[i]])
				if (sum (y>0) < 5){
					results <- append(results, list('zero'))
				}
				else if (rsem_table$length[[i]] < cutoff){
					print('transcript with low length')
					results <- append(results, list('zero'))
				} else {
					result <- fitdistr(y, 'negative binomial')
					results <- append(results, list(result))
				}
			} else if (model == 'lognormal') {
				if (sum (x[[i]] > 0) < 5) {
					results <- append(results, list('zero'))
				} else {
					y <- log(x[[i]] + 1)
					result <- glm(y~1, family=gaussian())
					results <- append(results, list(result))
				}
			} else {
				stop("does not match model")
			}
		},
		error = function(e)
		{
			print(e)
			print(i)
			results <<- append(results, list('zero'))
		})
	}
	names(results) <- names(x)
	results
}
#choose genes with only all nonzerotranscripts
choose_genes_to_perturb <- function(fits, fraction)
{
	transcript_names <- names(fits)
	zeros <- sapply(fits, function(x) typeof(x) != "list")
	zero_transcripts <- transcript_names[zeros]

	mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
	t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",  "external_gene_name"), mart = mart)
	t2g <- t2g[t2g$ensembl_transcript_id %in% names(fits),] 
	t2g_zero <- t2g[t2g$ensembl_transcript_id %in% zero_transcripts,]
	zero_genes <- unique(t2g_zero$external_gene_name)

	# select only genes without zero transcript
	t2g_nonzero <- t2g[!t2g$external_gene_name %in%  zero_genes, ]	
	nonzero_genes <- unique(t2g_nonzero$external_gene_name)
	num_nonzero_genes <- length(nonzero_genes)
	num_sample <- round(num_nonzero_genes * fraction)
	indices <- sample(1:num_nonzero_genes, num_sample, replace=FALSE)
	effect_sizes <- draw_effects_size(num_sample)
	perturbed_genes <- nonzero_genes[indices]

	#only perturb genes that are nonzero
	#transcript to gene table
	t2g_perturb <- t2g_nonzero[t2g_nonzero$external_gene_name%in% perturbed_genes, ]
	perturb_gene_table <- data.frame(perturbed_genes = perturbed_genes, effect_sizes = effect_sizes)
	t2g_perturb$perturbed_transcripts <- t2g_perturb$ensembl_transcript_id
	t2g_perturb$perturbed_genes<- t2g_perturb$external_gene_name
	t2g_perturb <- merge(t2g_perturb, perturb_gene_table, by=c('perturbed_genes'), left.all = TRUE) 
	#perturb_table required perturbed_transcripts and effect_sizes
	t2g_perturb
}


#select nonzero transcripts to perturb, return indices based on list of fits
choose_transcripts_to_perturb <- function(fits, fraction)
{
	all_transcripts <- names(fits)
	#select nonzero transcripts
	nonzeros <- sapply(fits, function(x) typeof(x) == "list")
	nonzero_transcripts <- all_transcripts[nonzeros]
	num_transcripts <- length(nonzero_transcripts)
	print(num_transcripts)
	num_sample <- round(num_transcripts * fraction)
	print(num_sample)
	#choose random indices to perturb
	indices <- sample(1:num_transcripts, num_sample, replace=FALSE)
	print(length(indices))
	perturbed_transcripts <- nonzero_transcripts[indices]
	indices <- match(perturbed_transcripts, all_transcripts)

	effect_sizes <- draw_effects_size(num_sample)

	#perturb_table required perturbed_transcripts and effect_sizes
	perturb_table <- data.frame(perturbed_transcripts, indices, effect_sizes)
	perturb_table <- perturb_table[order(indices),]

	perturb_table
}

#to change between log and nonlog model, exp versus do not exp the model. but both equate to lognormal effects with at least 2
#drawing effects size from normal
draw_effects_size <- function(nsample)
{
	#DESeq drew from log normal with mean =0
	#rlnorm(nsample, meanlog=0, sdlog=1.5)
	# I will do the same but now with truncated log normal
	#change to normal because in log counts
	x <- rtruncnorm(nsample, a=log(2), b=Inf, mean=0, sd=1)
	sign <- sample(c(-1,1), nsample, replace=TRUE)
	x <- x*sign
	x
}




simulate <- function(fits, ncells, perturb_table, model='negative binomial')
{
	#check if needs to perturb before simulating
	if(model == 'lognormal')
	{
		sims <- lapply(1:length(fits), function(i) {
			simulate_from_lognormal(fits[[i]], ncells, names(fits[i]), perturb_table)})
	}
	else if(model =='negative binomial')
	{
		sims <- lapply(1:length(fits), function(i) {
			simulate_from_nbinom(fits[[i]], ncells, names(fits[i]), perturb_table)})
	}
	else
	{
		stop('does not match model')
	}

	expect_equal(length(fits), length(sims))
	names(sims) <- names(fits)
	sims
}


#to change between log and nonlog model, change effectsize*mean versus effectsize + mean. do not need to retransform
simulate_from_lognormal<- function(fit, nsamples, transcript_id, perturb_table)
{
	#if not filtered out
	if(typeof(fit) == 'list') {
		mean <- fit$coefficients[[1]]
		sd <- sqrt(sum(fit$residual^2)/fit$df.residual)
		if(transcript_id %in% perturb_table$perturbed_transcripts) {
			index <- which(perturb_table$perturbed_transcripts == transcript_id)
			effect_size <- perturb_table$effect_sizes[index]
			samples <- rnorm(nsamples, mean=mean+effect_size, sd=sd)
			samples <- as.vector(exp(samples)-1)
			samples[samples<0] <- 0 
			return(samples)
		}
		else {
			samples <- rnorm(nsamples, mean=mean, sd=sd)
			samples <- as.vector(exp(samples) -1)
			samples[samples<0] <- 0
			return(samples)
		}
	}
	else{
		return(rep(0, nsamples))
	}
}

#to change between log and nonlog model, change effectsize*mean versus effectsize + mean. do not need to retransform
simulate_from_nbinom <- function(fit, nsamples, transcript_id, perturb_table)
{
	#if not filtered out
	if(typeof(fit) == 'list') {
		size <- fit$estimate[[1]]
		mu <- fit$estimate[[2]]
		if(transcript_id %in% perturb_table$perturbed_transcripts) {
			index <- which(perturb_table$perturbed_transcripts == transcript_id)
			effect_size <- perturb_table$effect_sizes[index]
			samples <- rnbinom(nsamples, mu = mu * exp(effect_size), size=size)
			return(samples)
		}
		else {
			samples <- rnbinom(nsamples, mu = mu , size= size)
			return(samples)
		}
	}
	else{
		return(rep(0, nsamples))
	}
}


#to reconvert the table into the format by cells instead of by transcripts
map_sims <- function(sims)
{
	ncells <- length(sims[[1]])
	tpms <- lapply(1:ncells, function(i) {sapply(sims, function(x) x[i])})
	tpms
}

# print table in the format that RSEM takes in to simulate reads
print_tables <- function(template_table, tpms, folder, type='tpms')
{
	ncells <- length(tpms)
	for(i in 1:ncells)
	{
		if(type == 'tpms')
			template_table$TPM <- tpms[[i]]
		else if (type == 'counts')
			template_table$count <- tpms[[i]]
		else
			stop ('cannot match type of reads')

		title <- toString(i)
		title <- paste(title, '.results', sep='')
		title <- paste(folder, title, sep='')
		write.table(template_table, file = title, quote=FALSE, row.names=FALSE, sep = '\t')
	
	}
}

#print table to external file
print_perturb <- function(perturb_table, path)
{
	write.table(perturb_table, file = path, quote=FALSE, row.names=FALSE, sep = '\t')
}
