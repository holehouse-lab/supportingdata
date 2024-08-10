#!/usr/bin/env python

'''
> make_seq_in_file.py FD.isf IDR.isf seq.in
'''

import sys
import argparse

parser = argparse.ArgumentParser(description='create seq.in file')
parser.add_argument('fd_seq', type=str, help='FD amino acid sequence')
parser.add_argument('idr_seqfile', type=str, help='single line IDR amino acid sequence file')
parser.add_argument('outfile', type=str, help='output file')
parser.add_argument('--idr-first', action='store_true', help='flag if IDR is first chain in PDB')
parser.add_argument('--idr-caps', action='store_true', help='cap IDR with ACE and NME')
parser.add_argument('--fd-caps', action='store_true', help='cap FD with ACE and NME')
args = parser.parse_args()

fd_seq = args.fd_seq
idr_seqfile = args.idr_seqfile
outfile = args.outfile

# Get IDR aa sequence from isf file
with open(idr_seqfile) as f:
    lines = [x.strip().split() for x in f]
    idr_seq = lines[0][1] 

# AA code conversion table
aa_code = {}
code1let = 'ACDEFGHIKLMNPQRSTVWY'
code3let = ['Ala','Cys','Asp','Glu','Phe','Gly','Hie','Ile','Lys','Leu',
            'Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr']
for i, aa in enumerate(list(code1let)):
    aa_code[aa] = code3let[i].upper()


if args.idr_first:
    # Iterate through IDR, adding residues to string
    if args.idr_caps:
        s = 'ACE\n' # Add IDR NTD cap
        s += f'{aa_code[idr_seq[0]]}\n'
        for idr_aa in list(idr_seq)[1:-1]:
            s += f'{aa_code[idr_aa]}\n'
        s += f'{aa_code[idr_seq[-1]]}\n'
        s += 'NME\n' # Add IDR CTD cap

    else:
        s = f'{aa_code[idr_seq[0]]}_N\n'
        for idr_aa in list(idr_seq)[1:-1]:
            s += f'{aa_code[idr_aa]}\n'
        s += f'{aa_code[idr_seq[-1]]}_C\n'

else:
    s = ''


# Iterate through FD, adding residues to string
if args.fd_caps: # Add ACE and NME caps to FD
    s += 'ACE\n'
    s += f'{aa_code[fd_seq[0]]}\n'
    for fd_aa in list(fd_seq)[1:-1]:
        s += f'{aa_code[fd_aa]}\n'
    s += f'{aa_code[fd_seq[-1]]}\n'
    s += 'NME\n'

else:
    s += f'{aa_code[fd_seq[0]]}_N\n'
    for fd_aa in list(fd_seq)[1:-1]:
        s += f'{aa_code[fd_aa]}\n'
    s += f'{aa_code[fd_seq[-1]]}_C\n'


if not args.idr_first:
    # Iterate through IDR, adding residues to string
    if args.idr_caps:
        s += 'ACE\n' # Add IDR NTD cap
        s += f'{aa_code[idr_seq[0]]}\n'
        for idr_aa in list(idr_seq)[1:-1]:
            s += f'{aa_code[idr_aa]}\n'
        s += f'{aa_code[idr_seq[-1]]}\n'
        s += 'NME\n' # Add IDR CTD cap

    else:
        s += f'{aa_code[idr_seq[0]]}_N\n'
        for idr_aa in list(idr_seq)[1:-1]:
            s += f'{aa_code[idr_aa]}\n'
        s += f'{aa_code[idr_seq[-1]]}_C\n'


# Get net charge
net_charge = 0
for i, fd_aa in enumerate(list(fd_seq)):
    if fd_aa in ['D', 'E']:
        net_charge -= 1
    elif fd_aa in ['R', 'K']:
        net_charge += 1
for idr_aa in list(idr_seq):
    if idr_aa in ['D', 'E']:
        net_charge -= 1
    elif idr_aa in ['R', 'K']:
        net_charge += 1

# Titrate ions to achieve neutral charge
while net_charge != 0:
    if net_charge < 0:
        s += 'NA+\n'
        net_charge += 1
    elif net_charge > 0:
        s += 'CL-\n'
        net_charge -= 1

s += 'END\n'

# Write to file
with open(outfile, 'w') as out:
    out.write(s)
