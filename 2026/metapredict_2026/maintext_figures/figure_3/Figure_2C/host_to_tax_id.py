import os
from pathlib import Path

# Inputs and the five per-host tables written here live in data/viral_proteomes/
# at the top of the repository (three levels above this directory).
DATA_DIR = Path(__file__).resolve().parents[3] / 'data' / 'viral_proteomes'

with open(DATA_DIR / 'virushostdb.tsv') as fh:
    lines=fh.read().split('\n')
fh.close()
headers=lines[0].split('\t')
eukaryotic=[]
archaea=[]
human=[]
bacteria=[]
ecoli=[]

# resolve column positions from the header once, rather than hard-coding them
virus_tax_id_idx = headers.index('virus tax id')
host_tax_id_idx = headers.index('host tax id')
host_lineage_idx = headers.index('host lineage')

for i in lines[1:]:
    fields = i.split('\t')
    # need enough columns to read the host lineage (its index is the largest of
    # the three we use, so this also guarantees the host/virus tax id columns are
    # present). guarding on >= 9 previously let short rows through and then raised
    # an IndexError on fields[9].
    if len(fields) > host_lineage_idx:
        host_kingdom = fields[host_lineage_idx].split(';')[0]
        if host_kingdom == 'Eukaryota':
            eukaryotic.append(fields[virus_tax_id_idx])
        elif host_kingdom == 'Archaea':
            archaea.append(fields[virus_tax_id_idx])
        elif host_kingdom == 'Bacteria':
            bacteria.append(fields[virus_tax_id_idx])
        if fields[host_tax_id_idx] == '9606':
            human.append(fields[virus_tax_id_idx])
        if fields[host_tax_id_idx] == '562':
            ecoli.append(fields[virus_tax_id_idx])


# get org id to prot id
org_id_to_prot_id={}
with open(DATA_DIR / 'org_id_to_proteome_id.tsv') as fh:
    lines=fh.read().split('\n')
fh.close()
for i in lines[1:]:
    if i !='':
        org_id_to_prot_id[i.split('\t')[1]]=i.split('\t')[0]

eukaryotic_proteomes=[]
bacterial_proteomes=[]
archaeal_proteomes=[]
human_proteomes=[]
ecoli_proteomes=[]


for i in eukaryotic:
    if i in org_id_to_prot_id:
        eukaryotic_proteomes.append(org_id_to_prot_id[i])

for i in archaea:
    if i in org_id_to_prot_id:
        archaeal_proteomes.append(org_id_to_prot_id[i])

for i in bacteria:
    if i in org_id_to_prot_id:
        bacterial_proteomes.append(org_id_to_prot_id[i])

for i in human:
    if i in org_id_to_prot_id:
        human_proteomes.append(org_id_to_prot_id[i])

for i in ecoli:
    if i in org_id_to_prot_id:
        ecoli_proteomes.append(org_id_to_prot_id[i])


# now get fractions of disorder. 
with open(DATA_DIR / 'viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv') as fh:
    lines=fh.read().split('\n')
org_id_to_disorder={}
for i in lines[1:]:
    if i !='':
        org_id_to_disorder[i.split('\t')[0]]={'per_protein':i.split('\t')[1], 'total_proteome':i.split('\t')[2]}
fh.close()

eukaryotic_disorder={}
bacterial_disorder={}
archaeal_disorder={}
human_disorder={}
ecoli_disorder={}


for i in eukaryotic_proteomes:
    if i in org_id_to_disorder:
        eukaryotic_disorder[i]=org_id_to_disorder[i]

for i in bacterial_proteomes:
    if i in org_id_to_disorder:
        bacterial_disorder[i]=org_id_to_disorder[i]

for i in archaeal_proteomes:
    if i in org_id_to_disorder:
        archaeal_disorder[i]=org_id_to_disorder[i]

for i in human_proteomes:
    if i in org_id_to_disorder:
        human_disorder[i]=org_id_to_disorder[i]

for i in ecoli_proteomes:
    if i in org_id_to_disorder:
        ecoli_disorder[i]=org_id_to_disorder[i]



# write out files
with open(DATA_DIR / 'eukaryotic_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in eukaryotic_disorder:
        fh.write(f"{i}\t{eukaryotic_disorder[i]['per_protein']}\t{eukaryotic_disorder[i]['total_proteome']}\n")

with open(DATA_DIR / 'bacterial_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in bacterial_disorder:
        fh.write(f"{i}\t{bacterial_disorder[i]['per_protein']}\t{bacterial_disorder[i]['total_proteome']}\n")

with open(DATA_DIR / 'archaeal_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in archaeal_disorder:
        fh.write(f"{i}\t{archaeal_disorder[i]['per_protein']}\t{archaeal_disorder[i]['total_proteome']}\n")

with open(DATA_DIR / 'human_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in human_disorder:
        fh.write(f"{i}\t{human_disorder[i]['per_protein']}\t{human_disorder[i]['total_proteome']}\n")

with open(DATA_DIR / 'ecoli_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in ecoli_disorder:
        fh.write(f"{i}\t{ecoli_disorder[i]['per_protein']}\t{ecoli_disorder[i]['total_proteome']}\n")

