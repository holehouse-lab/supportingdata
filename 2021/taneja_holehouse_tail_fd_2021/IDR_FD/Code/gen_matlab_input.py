import os
import sys 
import argparse
import numpy as np  
import pandas as pd
import math
from pathlib import Path
import Bio.PDB as bpdb
from Bio.PDB import PDBParser, PDBIO, Select  
from pathlib import Path
import glob
import random

def get_pdb_stem(pdb_file):

    if 'clean' not in pdb_file:
        out = (pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('.')])
    else:
        out = (pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('_clean')])

    #pdb2pqr has trouble if | is present in pdb_file
    if '|' in pdb_file:
        return(out[0:out.find('|')])
    else:
        return(out)


def get_atom_coords(pdb_file, chain_num, residue_type, custom_residue=-1):

    #residue_type:
    #if Nterm use min residue
    #if Cterm use max residue
    #if Custom use custom residue
    #if Random use random residue


    if residue_type == 'Custom':
        if custom_residue == -1:
            print("must pass in custom residue if nterm_or_cter not Nterm/Cterm")
            sys.exit()

    residue_list = []
    atom_num_list = []
    structure = PDBParser(QUIET=True).get_structure('temp', pdb_file)
    for model in structure:
        for chain in model:
            chain_name = (chain.get_full_id()[2])
            if chain_name == chain_num:
                for residue in chain:
                    residue_num = residue.get_full_id()[3][1]
                    residue_list.append(residue_num)

    min_residue = min(residue_list)
    max_residue = max(residue_list)

    num_residues = len(residue_list)


    if residue_type == 'N_terminal':
        relevant_residue = [min_residue]
    elif residue_type == 'C_terminal':
        relevant_residue = [max_residue]
    elif residue_type == 'Custom':
        if custom_residue in residue_list:
            relevant_residue = [custom_residue]
        else:
            print('%i not in PDB for %s' % (custom_residue, pdb_file))
            nearest_residue = min(residue_list, key=lambda x:abs(x-custom_residue))
            relevant_residue = [nearest_residue]
            print('Using residue %i instead' % nearest_residue)
    elif residue_type == 'Random':
        num_residues_to_sample = min(1 + max(0,(round(num_residues/100)-1)*1),4)
        relevant_residue = random.sample(residue_list, num_residues_to_sample)


    residue_coord = {}

    for curr_residue in relevant_residue:
        atom_num_list = []
        structure = PDBParser(QUIET=True).get_structure('temp', pdb_file)
        for model in structure:
            for chain in model:
                chain_name = (chain.get_full_id()[2])
                if chain_name == chain_num:
                    for residue in chain:
                        residue_num = residue.get_full_id()[3][1]
                        if residue_num == curr_residue:
                            for atom in residue:
                                atom_type = atom.get_name()
                                if atom_type == 'CA':
                                    atom_num = atom.get_serial_number()
                                    coords = list(residue[atom_type].get_vector())
                                    coords[0] = round(coords[0],2)
                                    coords[1] = round(coords[1],2)
                                    coords[2] = round(coords[2],2)
                                    residue_coord[curr_residue] = [coords[0],coords[1],coords[2]]

    if len(residue_coord.keys()) == 0:
        print('No residues selected')
        print(pdb_file)
        print(residue_type)
        print(custom_residue)
        print(residue_list)
        sys.exit()

    return(residue_coord)


apbs_input_dir = '../PDB_Structure_FD_candidates_within10_FD/APBS_Input/Ion_Concentration_10' 
apbs_output_dir = '../PDB_Structure_FD_candidates_within10_FD/APBS_Output/Ion_Concentration_10'
matlab_metadata_filename = '../Processed_Data/matlab_metadata.csv'

matlab_metadata = pd.read_csv(matlab_metadata_filename)

pdb_files = [file for file in glob.glob("%s/*.pdb" % apbs_input_dir)]
ply_files = [file for file in glob.glob("%s/*.ply" % apbs_output_dir)]
dx_files = [file for file in glob.glob("%s/*.dx" % apbs_output_dir)]



uniprotid_list = []
pdbid_list = []
chain_list = []
refpoint_list = []
xloc_list = []
yloc_list = []
zloc_list = []
dxloc_list = []
plyloc_list = []
res_list = []


for index, row in matlab_metadata.iterrows():

    print("On index %i" % index)

    uniprot_id = row['Uniprot_ID']
    pdb_id = row['PDB_ID']
    chain = row['chain']
    fd_start_real = row['fd_start_real']
    fd_start_pdb = row['fd_start_pdb']
    tail_type = row['tail']

    pdb_w_chain = pdb_id + '_' + chain

    pdb_file = '%s/%s_clean.pdb' % (apbs_input_dir, pdb_w_chain)
    dx_file = '%s/%s.dx' % (apbs_output_dir, pdb_w_chain)
    ply_file = '%s/%s.ply' % (apbs_output_dir, pdb_w_chain)

    atom_coords = get_atom_coords(pdb_file, chain, 'Custom', fd_start_pdb)

    for res in atom_coords:

        uniprotid_list.append(uniprot_id)
        pdbid_list.append(pdb_id)
        chain_list.append(chain)
        refpoint_list.append(tail_type)
        plyloc_list.append(ply_file)
        dxloc_list.append(dx_file)

        coords = atom_coords[res]
        res_list.append(res)
        xloc_list.append(coords[0])
        yloc_list.append(coords[1])
        zloc_list.append(coords[2])


    atom_coords = get_atom_coords(pdb_file, chain, 'Random', fd_start_pdb)
    
    for res in atom_coords:

        uniprotid_list.append(uniprot_id)
        pdbid_list.append(pdb_id)
        chain_list.append(chain)
        refpoint_list.append('Random')
        plyloc_list.append(ply_file)
        dxloc_list.append(dx_file)

        coords = atom_coords[res]
        res_list.append(res)
        xloc_list.append(coords[0])
        yloc_list.append(coords[1])
        zloc_list.append(coords[2])



    out = pd.DataFrame({'Uniprot_ID': uniprotid_list, 'PDB_ID': pdbid_list, 'chain': chain_list,
                        'reference_point': refpoint_list, 'residue': res_list, 'x_loc': xloc_list, 'y_loc': yloc_list, 'z_loc': zloc_list,
                        'ply_file_loc': plyloc_list, 'dx_file_loc': dxloc_list})

    
out.to_csv('../Processed_Data/phi_patchiness_input.csv', index=False)
print(out)


