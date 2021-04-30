#pdb2pqr --> pqr + in
#apbs --> dx
#ETSurf --> ply 

import sys
import os 
import subprocess
import fileinput
from pathlib import Path



#insert executable path to pdb2pqr, apbs, and EDTsurf

pqr_path = '/Users/ishan/Holehouse/pdb2pqr-osx-bin64-2.1.0/pdb2pqr'
apbs_path = '/Users/ishan/Holehouse/apbs/build/bin/apbs'
edtsurf_path = '/Users/ishan/Holehouse/EDTSurf_mac'


def get_pdb_stem(pdb_file):

	if 'clean' not in pdb_file:
		return(pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('.')])
	else:
		return(pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('_clean')])


def clean_pdb_file_for_pqr(pdb_file, apbs_input_dir):

	pdb_clean_file = apbs_input_dir + '/' + get_pdb_stem(pdb_file) + '_clean.pdb'

	with open(pdb_file, "r") as inputFile,open(pdb_clean_file,"w") as outFile:
		for line in inputFile:
			if line.startswith("ANISOU") == False and line.startswith("SEQADV") == False and line.startswith("MDLTYP") == False and line.startswith("HETATM") == False and line.startswith('CRYST1') == False and line.startswith('DBREF1') == False and line.startswith('DBREF2') == False:
				outFile.write(line)

	return(pdb_clean_file)


def convert_pdb_to_pqr_amber(pdb_clean_file, apbs_input_dir):

	pdb_stem = get_pdb_stem(pdb_clean_file)
	pqr_file = '%s/%s.pqr' % (apbs_input_dir, pdb_stem)

	Path(apbs_input_dir).mkdir(parents=True, exist_ok=True)

	pqr_cmd = pqr_path + ' --apbs-input --ff=AMBER ' + pdb_clean_file + ' ' + pqr_file
	subprocess.call(pqr_cmd, shell=True)
	subprocess.call('rm ' + apbs_input_dir + '/' + pdb_stem + '-input.p', shell=True)


def convert_pdb_to_pqr(pdb_clean_file, apbs_input_dir):

	pdb_stem = get_pdb_stem(pdb_clean_file)
	pqr_file = '%s/%s.pqr' % (apbs_input_dir, pdb_stem)

	Path(apbs_input_dir).mkdir(parents=True, exist_ok=True)

	pqr_cmd = pqr_path + ' --apbs-input --ff=PARSE ' + pdb_clean_file + ' ' + pqr_file
	subprocess.call(pqr_cmd, shell=True)
	subprocess.call('rm ' + apbs_input_dir + '/' + pdb_stem + '-input.p', shell=True)




def update_apbs_in(pdb_file, apbs_input_dir, temp, ion_conc):

	pdb_stem = get_pdb_stem(pdb_file)
	apbs_input_file = '%s/%s.in' % (apbs_input_dir, pdb_stem)

	pqr_file = '%s/%s.pqr' % (apbs_input_dir, pdb_stem)

	apbs_output_dir = apbs_input_dir.replace('APBS_Input', 'APBS_Output')
	dx_file_stem = '%s/%s' % (apbs_output_dir, pdb_stem)


	if os.path.exists(apbs_input_file):
		with fileinput.FileInput(apbs_input_file, inplace=True) as file:
			for line in file:
				if 'mol pqr' in line:
					line = "    mol pqr %s\n" % pqr_file
				elif 'write pot dx' in line:
					line = "    write pot dx %s\n" % dx_file_stem
				elif 'temp' in line:
					line = "    temp %i" % temp
					if ion_conc is not None:
						line = line + '\n' + "    ion charge +1 conc %.1g radius 2.0\n    ion charge -1 conc %.1g radius 1.8\n" % (ion_conc, ion_conc)

				sys.stdout.write(line)


def gen_pqr_apbs_input(pdb_file, apbs_input_dir, temp, ion_conc=None):

	pdb_clean_file = clean_pdb_file_for_pqr(pdb_file, apbs_input_dir)
	convert_pdb_to_pqr(pdb_clean_file, apbs_input_dir)
	update_apbs_in(pdb_file, apbs_input_dir, temp, ion_conc)



def gen_pqr_apbs_input_amber(pdb_file, apbs_input_dir):

	pdb_clean_file = clean_pdb_file_for_pqr(pdb_file, apbs_input_dir)
	convert_pdb_to_pqr_amber(pdb_clean_file, apbs_input_dir)
	update_apbs_in(pdb_file, apbs_input_dir)

def gen_apbs_dx(apbs_in_file):

	apbs_cmd = apbs_path + ' ' + apbs_in_file
	subprocess.call(apbs_cmd, shell=True)



def gen_edtsurf_output(pdb_file, apbs_input_dir):

	#if 'clean' not in pdb_file:
	#	print("not a clean PDB file")
	#	sys.exit()

	pdb_stem = get_pdb_stem(pdb_file)

	apbs_output_dir = apbs_input_dir.replace('APBS_Input', 'APBS_Output') + '/' + pdb_stem

	edtsurf_cmd = edtsurf_path + ' -i ' + pdb_file + ' -o ' + apbs_output_dir
	subprocess.call(edtsurf_cmd, shell=True)




