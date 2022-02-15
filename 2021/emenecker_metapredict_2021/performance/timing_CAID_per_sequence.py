from metapredict import meta
import protfasta
import random
import time
import numpy as np

## This script calculates the per-sequence time for each of the sequences in the DisProt dataset. We ultimately do not 
## use this data in the paper, but include this anyway because it shows the speed with which metapredict can compute
## disorder and may be useful for someone else
##

# read in the CAID sequences using protfasta

# Note - if you do not have protfasta it can be installed using
#    pip install protfasta
#
# Documentation available here:
# https://protfasta.readthedocs.io/en/latest/

LEGACY=True

print('Reading in CAID sequences...')
disprot_data = protfasta.read_fasta('../data/CAID/CAID_DisProt_Dataset.fasta')

# caid_id_order is just the order of IDs as reported in https://codeocean.com/capsule/2223745/tree/v1/dataset_stats/cpu_times.csv
with open('data/caid_id_order.csv','r') as fh:
    content = fh.readlines()

ordered_ids = [i.strip() for i in content]

number_of_iterations = 3

per_protein_time = {}
print('Running loops')

c= 0
for idx in ordered_ids:
    c=c+1
    
    per_protein_time[idx] = []
    s = disprot_data[idx]
    print('On %i of %i' % (c, len(ordered_ids)))

    for i in range(number_of_iterations):        

        start = time.time()
        cur_dis = meta.predict_disorder(s, legacy = LEGACY)
        end = time.time()
        loop_time = end - start
        per_protein_time[idx].append(loop_time)

    per_protein_time[idx] = np.mean(per_protein_time[idx])


print('DONE')

if LEGACY is True:
    with open('per_protein_times.csv', 'w') as fh:
        for idx in ordered_ids:
            fh.write(f'{idx}, {per_protein_time[idx]}\n')

else:
    with open('per_protein_times_metapredict_v2.csv', 'w') as fh:
        for idx in ordered_ids:
            fh.write(f'{idx}, {per_protein_time[idx]}\n')
