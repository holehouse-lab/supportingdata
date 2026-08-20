#!/usr/bin/env python

'''
modify_complex_pdb.py

Provide a full PDB file of a FD-IDR complex and output a cleaned version with superfluous residues removed.
Particularly useful for converting a AF2-predicted complex for integration into simulation pipeline.

Usage:

> python modify_complex_pdb.py <arg1> <arg2> <arg3>


'''

import argparse
import os
import numpy as np
from housetools.parameters import AADICT_3to1

## --------------------- Functions --------------------- ##
# Add to this set of keywords as needed
KEEP_KEYWORDS = ['MODEL', 'ATOM', 'TER', 'ENDMDL', 'END']

def read_input_pdb(pdbfile):
    # TODO: modify this function to work in cases where there are >1000 residues
    # - so split by position in string rather than by whitespace delimiters
    if os.path.exists(pdbfile):
        with open(pdbfile) as f:
            lines = [x.strip().split() for x in f]

        # First pass: remove all lines that do not start with relevant keywords
        clean_lines = [line for line in lines if line[0] in KEEP_KEYWORDS]

        # Return this list of lines to be filtered further
        return clean_lines

    else:
        raise Exception(f'Provided PDB file does not exist: {pdbfile}')

def get_atom_line_info(pdb_lines):
    atom_line_idxs = np.array([i for i in range(len(pdb_lines)) if pdb_lines[i][0] == 'ATOM'], dtype=int)
    atom_lines = [pdb_lines[i] for i in atom_line_idxs]

    residues = []
    chain_ids = []
    res_idxs = []
    prev_res_idx = -1
    for l in atom_lines:
        res_name = l[3]
        chain_id = l[4]
        cur_res_idx = l[5]

        if cur_res_idx != prev_res_idx:
            try:
                residues.append(AADICT_3to1[res_name])
            except:
                if res_name == 'ACE' or res_name == 'NME':
                    residues.append('=')
            chain_ids.append(chain_id)
            res_idxs.append(cur_res_idx)

        prev_res_idx = cur_res_idx

    # Sort residues by chain
    sequences = []
    chain_res_idxs = []
    cur_sequence = []
    cur_chain_res_idx = []
    prev_chain_id = ''
    for aa, chain, res_idx in zip(residues, chain_ids, res_idxs):
        if chain != prev_chain_id:
            sequences.append("".join(cur_sequence))
            chain_res_idxs.append(cur_chain_res_idx)
            cur_sequence = [aa]
            cur_chain_res_idx = [int(res_idx)]
        else:
            cur_sequence.append(aa)
            cur_chain_res_idx.append(int(res_idx))

        prev_chain_id = chain

    sequences.append("".join(cur_sequence))
    chain_res_idxs.append(cur_chain_res_idx)

    sequences = sequences[1:]
    chain_res_idxs = chain_res_idxs[1:]
    chains = np.unique(chain_ids)

    return chains, sequences, chain_res_idxs

def edit_TER_line(preceding_pdb_line):
    ter_atom_idx = str(int(preceding_pdb_line[1]) + 1)
    ter_res_name = preceding_pdb_line[3]
    ter_res_chain = preceding_pdb_line[4]
    ter_res_idx = preceding_pdb_line[5]

    return ['TER', ter_atom_idx, ter_res_name, ter_res_chain, ter_res_idx]



## --------------------- MAIN --------------------- ##

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pdb', required=True, help='Input PDB file')
parser.add_argument('--out', default=None, help='Output PDB file (default is "out_<inputPDBname>")')
parser.add_argument('--clip-idx', nargs=2, type=int, default=None, help='Residue index range to keep from IDR chain (inclusive)')
parser.add_argument('--clip-str', type=str, default=None, help='NOT IMPLEMENTED YET')
parser.add_argument('--idr-first', action='store_true', default=False)
args = parser.parse_args()

pdb_lines = read_input_pdb(args.pdb)
chains, sequences, res_idxs = get_atom_line_info(pdb_lines)

if args.idr_first:
    idr_chain = chains[0]
    idr_seq = sequences[0]
else:
    idr_chain = chains[1]
    idr_seq = sequences[1]

if len(chains) > 2:
    print(f'PDB with >2 chains provided. Ignoring chains {chains[2:]}...')

# Clean up unwanted residues in PDB
if args.clip_str is not None and args.clip_idx is not None:
    raise Exception('Cannot use --clip-idx and --clip-str together')

elif args.clip_idx is not None:
    clipped_pdb_lines = []

    for line in pdb_lines:
        # Add non-atom lines, with special case being the TER line of the clipped chain
        if line[0] != 'ATOM':
            if line[0] == 'TER':
                clipped_pdb_lines.append(edit_TER_line(clipped_pdb_lines[-1]))

            else:
                clipped_pdb_lines.append(line)

        # Add ATOM line if not IDR chain, or if IDR chain and within specified residue indices
        else:
            if line[4] != idr_chain:
                clipped_pdb_lines.append(line)
            else:
                if int(line[5]) >= args.clip_idx[0] and int(line[5]) <= args.clip_idx[1]:
                    clipped_pdb_lines.append(line)

elif args.clip_str is not None:
    clipped_pdb_lines = []
    # TODO implement (maybe find corresponding res_idxs first and use clip_idx code ^)

else:
    clipped_pdb_lines = pdb_lines


# Create PDB string
out_str = ''
for l in clipped_pdb_lines:
    if l[0] == 'ENDMDL' or l[0] == 'END':
        s = f'{l[0]}\n'
    elif l[0] == 'MODEL':
        s = f'MODEL {str.rjust(l[1], 5)}\n'
    elif l[0] == 'TER':
        s = f'TER   {str.rjust(l[1],5)}      {l[2]} {l[3]}{str.rjust(l[4],4)}\n'

    elif l[0] == 'ATOM':
        s = f'ATOM  {str.rjust(l[1],5)}  {str.ljust(l[2],4)}{l[3]} {l[4]}{str.rjust(l[5],4)}    '

        if len(l) > 11:
            s += f'{format(float(l[6]),"8.3f")}{format(float(l[7]),"8.3f")}{format(float(l[8]),"8.3f")}{format(float(l[9]),"6.2f")}{format(float(l[10]),"6.2f")}           {l[11]}\n'
        else: 
            s += f'{format(float(l[6]),"8.3f")}{format(float(l[7]),"8.3f")}{format(float(l[8]),"8.3f")}{format(float(l[9]),"6.2f")}{format(float(l[10]),"6.2f")} \n'

    else:
        raise Exception(f'Unrecognized keyword: {l[0]}')
    
    out_str += s

# Write to file
if args.out is not None:
    out_file_name = args.out
else:
    path = args.pdb.split('/')
    out_file_name = "/".join(path[:-1]) + "/out_" + path[-1]

with open(out_file_name, 'w') as out_f:
    out_f.write(out_str)


