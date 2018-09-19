import os
import multiprocessing as mp
import itertools
import sys

dir1 = sys.argv[1]
dir2 = sys.argv[2]
batchfile = sys.argv[3]

flnames = os.listdir(dir1)
flnames = [x for x in flnames if x.endswith('.fq')]
flnames = [int(x.split('.')[0]) for x in flnames]
flnames = sorted(list(set(flnames)))

flnames2 = os.listdir(dir2)
flnames2 = [x for x in flnames2 if x.endswith('.fq')]
flnames2 = [int(x.split('.')[0]) for x in flnames2]
flnames2 = sorted(list(set(flnames2)))

print(len(flnames))
print(len(flnames2))
f = open(batchfile, 'w')
for flname in flnames:
    flname1 = dir1 + str(flname) + '.results_1.fq'
    flname2 = dir1 + str(flname) + '.results_2.fq'
    f.write(str(flname) + '\t' + flname1 + '\t' + flname2 + '\n')
for flname in flnames2:
    flname1 = dir2 + str(flname) + '.results_1.fq'
    flname2 = dir2 + str(flname) + '.results_2.fq'
    f.write(str(flname) + '\t' + flname1 + '\t' + flname2 + '\n')

