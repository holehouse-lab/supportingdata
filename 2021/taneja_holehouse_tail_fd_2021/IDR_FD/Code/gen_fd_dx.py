import sys
import os
import gen_surface_charge_distribution
from os.path import isfile, join
from multiprocessing import Process
from concurrent.futures import ThreadPoolExecutor
import numpy as np

ion_conc = .10
ion_conc_str = '10'

apbs_input_dir = '../PDB_Structure_FD_candidates_within10/APBS_Input/Ion_Concentration_%s' % (ion_conc_str) 
apbs_in_files = sorted([apbs_input_dir + '/' + f for f in os.listdir(apbs_input_dir) if isfile(join(apbs_input_dir, f)) and '.in' in f])

def create_dx_parallel(apbs_in_files):
	for i in range(0,len(apbs_in_files)):
		apbs_in = apbs_in_files[i]
		print(i)
		print(apbs_in)
		gen_surface_charge_distribution.gen_apbs_dx(apbs_in)

def run_parallel(apbs_in_files_chunk):
    
    processes = []
    for i in range(0,len(apbs_in_files_chunk)):
        print('thread' + str(i))
        proc = Process(target=create_dx_parallel, args = (apbs_in_files_chunk[i],))
        proc.start()
        processes.append(proc)
    for proc in processes:
        proc.join()


def run_parallel_wrapper():
	apbs_in_files_chunks_7 =  list(np.array_split(np.array(apbs_in_files),7))
	run_parallel(apbs_in_files_chunks_7)

if __name__ == '__main__':
	run_parallel_wrapper()
