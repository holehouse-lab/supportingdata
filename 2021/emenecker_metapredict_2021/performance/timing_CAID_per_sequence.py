from metapredict import meta
import protfasta
import random
import time
import numpy as np

# read in the CAID sequences using protfasta

# Note - if you do not have protfasta it can be installed using
#    pip install protfasta
#
# Documentation available here:
# https://protfasta.readthedocs.io/en/latest/

print('Reading in CAID sequences...')
disprot_data = protfasta.read_fasta('../data/CAID/CAID_DisProt_Dataset.fasta')

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
        cur_dis = meta.predict_disorder(s)
        end = time.time()
        loop_time = end - start
        per_protein_time[idx].append(loop_time)

    per_protein_time[idx] = np.mean(per_protein_time[idx])


print('DONE')
with open('per_protein_times.csv', 'w') as fh:
    for idx in ordered_ids:
        fh.write(f'{idx}, {per_protein_time[idx]}\n')
