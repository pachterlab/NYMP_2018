# Ran Jan 16 to simulate reads with RSEM from the monocle model
# The monocle model is a tobit model on the cluster1 log tpm
# perturbations made on 20% of transcripts with at least 2 fold change, simulated from a log normal effect size distribution
# also used to simulate nonperturb

import os
import glob
import numpy
import sys

RSEM_path="/home/lynnyi/RSEM-1.3.0"
RSEM_command= RSEM_path + "/rsem-simulate-reads"
RSEM_ref_path= RSEM_path + "/ref/human"
RSEM_model= RSEM_path + "/exp/Trapnell_cluster1.stat/Trapnell_cluster1.model"
theta = "0.2"
i = sys.argv[1]
o = sys.argv[2]

reads_logmean = 14.42
reads_logsd = 0.3336 

with open('/home/lynnyi/dirichlet/log.txt', 'w') as logfile:
    reads = int(numpy.random.lognormal(mean=reads_logmean, sigma = reads_logsd))
    #reads = 2000000 
    cmd = RSEM_command + " " + RSEM_ref_path + " " + RSEM_model + " " + i + " " + theta + " " + str(reads) + " " + o + "  --seed 0 "
    logfile.write(cmd + '\n')
    print(cmd)
    os.system(cmd)

