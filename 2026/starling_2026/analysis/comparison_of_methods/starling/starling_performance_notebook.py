#
# Starling Performance Benchmarking Notebook; this notebook generates a few files:
#
# - results.pkl, a pickle file containing radius of gyration and time per conformation 
#    for each sequence, as well as the sequences themselves
#
# - diagnostics.csv, a csv file with sequence name, sequence, and time per conformation
#
# - starling_performance_time_per_conf_bs_<batch_size>_<compiled/not_compiled>.csv, a 
#   csv file with two columns: sequence length and time per conformation
#
import starling
from starling import generate
import protfasta
import time
import numpy as np
import pickle

# flag to set if we use the compile option
USE_COMPILE = False

if USE_COMPILE:
    starling.set_compilation_options(enabled=True)

batch_size = 100
n_conformations = 400
seq = protfasta.read_fasta('all_comparison_seqs.fasta')

rg_values = {}
ensemble_time = {}

for k in seq:

    start = time.time()
    e = generate(seq[k], conformations=n_conformations,return_single_ensemble=True, batch_size=batch_size)    
    end = time.time()
    ensemble_time[k] = (end-start)/n_conformations
    rg_values[k] = e.radius_of_gyration(return_mean=True)

# save out for diagnostics
diagnostics_data = {"rg_values": rg_values, "ensemble_time": ensemble_time, "sequence": seq}    
with open("results.pkl", "wb") as f:
    pickle.dump(diagnostics_data, f)

with open('diagnostics.csv','w') as fh:
    for k in ensemble_time:
        fh.write(f"{k}, {seq[k]}, {ensemble_time[k]}\n")
    
seq_len = []
ensemble_time_list = []
for k in seq:
    seq_len.append(len(seq[k]))
    ensemble_time_list.append(ensemble_time[k])


paired = list(zip(seq_len, ensemble_time_list))
paired_sorted = sorted(paired, key=lambda x: x[0])

# extract back
seq_len_sorted, time_per_conf_sorted = zip(*paired_sorted)

if USE_COMPILE:
    compile_string='compiled'
else:
    compile_string = 'not_compiled'
np.savetxt(f'starling_performance_time_per_conf_bs_{batch_size}_{compile_string}.csv', np.array([seq_len_sorted, time_per_conf_sorted]))
