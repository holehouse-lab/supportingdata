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

number_of_iterations = 3

loops = []
print('Running loops')
for i in range(0, number_of_iterations):
    print('On loop %i or %i' % (i+1, number_of_iterations))

    # start timer
    start = time.time()

    # cycle through each sequence and predict disorder (note we're not saving it!)
    for idx in disprot_data:
        s = disprot_data[idx]
        cur_dis = meta.predict_disorder(s)

    # end timer
    end = time.time()
    loop_time = end - start

    # append results to loops
    loops.append(loop_time)
    # print result
    print("==> metapredict takes %i seconds to go through the CAID DisProt Dataset"%(loop_time))

number_of_residues = 0
for idx in disprot_data:
    s = disprot_data[idx]
    number_of_residues = number_of_residues + len(s)


# print mean and standard error
mean_time = np.mean(loops)
std_error = np.std(loops)/np.sqrt(number_of_iterations)

print('Mean time = %i seconds +/- %i (stderr)'%(mean_time, std_error))
print('Average of %i residues per second' % (number_of_residues/mean_time))

