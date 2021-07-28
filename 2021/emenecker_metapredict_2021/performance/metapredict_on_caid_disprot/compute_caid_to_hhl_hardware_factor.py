## This code takes per-sequence disorder prediction times computed on our 
## local holehouse lab (HHL) hardware and calculates the ratio against per
## sequence calcualtions calculated on the CAID hardware, and then prints
## out the average. Predictions done using AUCPreD because it was much easier
## to install than any of the other top-performing predictors!
##
##


import numpy as np

# read in caid times 
with open('../data/caid_AUCPreD_times.csv') as fh:
    content = fh.readlines()
caid_times = {}
for line in content:
    sl = line.strip().split(',')
    caid_times[sl[0]] = float(sl[1])

# read in hhl times
with open('../data/aucpred_vs_metapredict_same_hardware_comparison.tsv') as fh:
    content = fh.readlines()

aucpred_hhl = {}
metapredict_hhl = {}

# [1:] so we skip the header line
for line in content[1:]:
    sl = line.strip().split()
    aucpred_hhl[sl[0]] = float(sl[5])
    metapredict_hhl[sl[0]] = float(sl[6])


ratios = []
for i in metapredict_hhl:
    if i in caid_times:
        ratios.append(aucpred_hhl[i]/caid_times[i])

print(f'Local HHL vs. CAID hardware conversion factor: {np.mean(ratios)}')
