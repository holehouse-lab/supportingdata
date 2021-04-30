import sys
import os
import gen_surface_charge_distribution
from os.path import isfile, join

ion_conc = .10
ion_conc_str = '10'

apbs_input_dir = '../PDB_Structure_FD_candidates_within10/APBS_Input/Ion_Concentration_%s' % (ion_conc_str) 
pdb_input_dir = apbs_input_dir
pdb_files = sorted([pdb_input_dir + '/' + f for f in os.listdir(pdb_input_dir) if isfile(join(pdb_input_dir, f)) and '.pdb' in f])

print(len(pdb_files))
for i in range(0,len(pdb_files)):
	pdb = pdb_files[i]
	print(i)
	print(pdb)
	gen_surface_charge_distribution.gen_edtsurf_output(pdb, apbs_input_dir)



