import os

with open(f'virushostdb.tsv') as fh:
    lines=fh.read().split('\n')
fh.close()
headers=lines[0].split('\t')
eukaryotic=[]
archaea=[]
human=[]
bacteria=[]
ecoli=[]

for i in lines[1:]:
    if len(i.split('\t'))>=9:
        if i.split('\t')[9].split(';')[0]=='Eukaryota':
            eukaryotic.append(i.split('\t')[0])
        elif i.split('\t')[9].split(';')[0]=='Archaea':
            archaea.append(i.split('\t')[0])
        elif i.split('\t')[9].split(';')[0]=='Bacteria':
            bacteria.append(i.split('\t')[0])
        if i.split('\t')[lines[0].split('\t').index('host tax id')] == '9606':
            human.append(i.split('\t')[0])
        if i.split('\t')[lines[0].split('\t').index('host tax id')] == '562':
            ecoli.append(i.split('\t')[0])


# get org id to prot id
org_id_to_prot_id={}
with open('org_id_to_proteome_id.tsv') as fh:
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
with open('viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv') as fh:
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
with open('eukaryotic_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in eukaryotic_disorder:
        fh.write(f"{i}\t{eukaryotic_disorder[i]['per_protein']}\t{eukaryotic_disorder[i]['total_proteome']}\n")

with open('bacterial_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in bacterial_disorder:
        fh.write(f"{i}\t{bacterial_disorder[i]['per_protein']}\t{bacterial_disorder[i]['total_proteome']}\n")

with open('archaeal_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in archaeal_disorder:
        fh.write(f"{i}\t{archaeal_disorder[i]['per_protein']}\t{archaeal_disorder[i]['total_proteome']}\n")

with open('human_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in human_disorder:
        fh.write(f"{i}\t{human_disorder[i]['per_protein']}\t{human_disorder[i]['total_proteome']}\n")

with open('ecoli_infecting_viruses.tsv', 'w') as fh:
    fh.write('proteome_id\tper_protein\ttotal_proteome\n')
    for i in ecoli_disorder:
        fh.write(f"{i}\t{ecoli_disorder[i]['per_protein']}\t{ecoli_disorder[i]['total_proteome']}\n")

