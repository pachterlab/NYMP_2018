source('simulation_framework.R')

args <- commandArgs(TRUE)
np_dir <- args[[1]]
p_dir <- args[[2]]
frac_perturb <- as.numeric(args[[3]])
model <- args[[4]]
sim_genes <- as.logical(args[[5]])

print('matching model')
if (model == 'negative binomial') {
	reads <- readRDS('~/dirichlet/hsmm_c1counts.rds')
	type <- 'counts'
} else if (model == 'lognormal') {
	reads <- readRDS('~/dirichlet/hsmm_c1tpms.rds')
	type <- 'tpms'
} else {
	stop("model does not match")
}

print('prepping template table')
template_path <- '/home/lynnyi/dirichlet/sims/test/RSEM/template.results'
template_table <- read.table(template_path, header=TRUE, sep='\t')

print('fitting')
fits <- fit(reads, model = model, rsem_table = template_table)
print(length(fits))
expect_equal(as.vector(template_table$transcript_id),names(fits))

print('perturbing table')
if (sim_genes) {
	perturb_table <- choose_genes_to_perturb(fits, frac_perturb)
} else { 
	perturb_table <- choose_transcripts_to_perturb(fits, frac_perturb)
}

print(dim(perturb_table))
print_perturb(perturb_table, file.path(p_dir, '/perturb.tsv'))

stopifnot(FALSE)

print('simulating')
sims_p <- simulate(fits, 105, perturb_table, model=model)
print(length(sims_p))

print('transposing table')
tpms_p <- map_sims(sims_p)
print(length(tpms_p))

print('printing sims')
print_tables(template_table, tpms_p, p_dir, type=type)

print('simulating null state')
sims_np<- simulate(fits, 105, NULL, model=model)

('transposing null table')
tpms_np<- map_sims(sims_np)

print('printing null table')
print_tables(template_table, tpms_np, np_dir, type=type) 

saveRDS(tpms_np, paste0(np_dir, '/tpms_np.rds'))
saveRDS(tpms_p, paste0(p_dir, '/tpms_p.rds'))

print('done')

