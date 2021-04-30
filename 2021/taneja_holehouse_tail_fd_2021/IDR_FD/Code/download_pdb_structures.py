#get PDB


import pandas as pd
import Bio
from Bio.PDB import PDBList, PDBParser, PDBIO, Select
import sys
from multiprocessing import Process
from concurrent.futures import ThreadPoolExecutor
import numpy as np

fd_data = pd.read_csv('../Processed_Data/fd_candidate.csv')
pdb_list = list(fd_data['PDB_ID'])
save_dir = '../PDB_Structure_FD_candidates_within10/PDB_Files'

'''pdbl = PDBList()
for i in pdb_list:
	pdbl.retrieve_pdb_file(i,file_format = 'pdb', pdir='../PDB_Structure_FD_candidates_within10/PDB_Files')
	print(i)
	io = PDBIO()
	pdb = PDBParser().get_structure(i, "%s/pdb%s.ent" % (save_dir,i.lower()))
	for chain in sorted(pdb.get_chains()):
		io.set_structure(chain)
		io.save(save_dir + '/' + pdb.get_id().lower() + "_" + chain.get_id() + ".pdb")
		break'''


def get_pdb_parallel(pdb_list):
	pdbl = PDBList()
	for i in pdb_list:
		pdbl.retrieve_pdb_file(i,file_format = 'pdb', pdir='../PDB_Structure_FD_candidates_within10/PDB_Files')
		print(i)
		io = PDBIO()
		pdb = PDBParser().get_structure(i, "%s/pdb%s.ent" % (save_dir,i.lower()))
		for chain in sorted(pdb.get_chains()):
			io.set_structure(chain)
			io.save(save_dir + '/' + pdb.get_id().lower() + "_" + chain.get_id() + ".pdb")

def run_parallel_wrapper(pdb_list_chunks):
    processes = []
    for i in range(0,len(pdb_list_chunks)):
        print('thread' + str(i))
        proc = Process(target=get_pdb_parallel, args = (pdb_list_chunks[i],))
        proc.start()
        processes.append(proc)
    for proc in processes:
        proc.join()





if __name__ == '__main__':
	pdb_list_chunks_7 = list(np.array_split(np.array(pdb_list),8))
	run_parallel_wrapper(pdb_list_chunks_7)




