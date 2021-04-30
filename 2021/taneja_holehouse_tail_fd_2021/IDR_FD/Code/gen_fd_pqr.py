import sys
import os
import gen_surface_charge_distribution
from os.path import isfile, join
import pandas as pd

temp = 310
ion_conc = .10
ion_conc_str = '10'

apbs_input_dir = '../PDB_Structure_FD_candidates_within10/APBS_Input/Ion_Concentration_%s' % (ion_conc_str) 
pdb_input_dir = '../PDB_Structure_FD_candidates_within10/PDB_Files'

matlab_fd_data = pd.read_csv('../Processed_Data/matlab_metadata.csv')
matlab_fd_data['PDB_ID_chain']=matlab_fd_data['PDB_ID']+'_'+matlab_fd_data['chain']
pdb_chain_list = sorted(list(set(list(matlab_fd_data['PDB_ID_chain']))))

pdb_files = ([pdb_input_dir + '/' + f + '.pdb' for f in pdb_chain_list])


for i in range(0,len(pdb_files)):
	pdb = pdb_files[i]
	print(i)
	print(pdb)
	gen_surface_charge_distribution.gen_pqr_apbs_input(pdb, apbs_input_dir, temp, ion_conc)


