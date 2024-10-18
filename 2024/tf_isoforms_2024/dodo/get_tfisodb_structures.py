# code for getting af2 structures from the internets
import requests
import io
from getSequence import getseq
import os
from dodo.build import pdb_from_pdb
from dodo.pdb_tools import PDBParser
from dodo.find_idrs_fds_loops import get_fds_idrs_from_metapredict

# function to download stuff from tfiso
def get_tfiso_pdb(protein_name, save_to=None):
    '''
    function to get the pdb from the tfisodb

    protein_name : string
        the protein name as a string. Needs to match the tfisodb name. 
    '''

    #generate the URL
    url=f'http://tfisodb.org/data/alphafold/{protein_name}/ranked_0.pdb'

    # download the af2 structure
    response=requests.get(url)
    in_memory_file = io.StringIO(response.text)
    lines=[]
    for line in in_memory_file:
        lines.append(line.strip())
    if save_to==None:
        return lines
    else:
        with open(save_to, 'w') as fh:
            for line in lines:
                fh.write(line+'\n')
        fh.close()
    
# download the structures
download_structures=False
if download_structures:
    get_tfiso_pdb('DLX4-1', save_to='dlx4_1.pdb')
    get_tfiso_pdb('DLX4-2', save_to='dlx4_2.pdb')
    get_tfiso_pdb('PKNOX1-1', save_to='pknox1_1.pdb')
    get_tfiso_pdb('PKNOX1-3', save_to='pknox1_3.pdb')

# set sequences so we can grab length. 
dlx4_1='MTSLPCPLPGRDASKAVFPDLAPVPSVAAAYPLGLSPTTAASPNLSYSRPYGHLLSYPYTEPANPGDSYLSCQQPAALSQPLCGPAEHPQELEADSEKPRLSPEPSERRPQAPAKKLRKPRTIYSSLQLQHLNQRFQHTQYLALPERAQLAAQLGLTQTQVKIWFQNKRSKYKKLLKQNSGGQEGDFPGRTFSVSPCSPPLPSLWDLPKAGTLPTSGYGNSFGAWYQHHSSDVLASPQMM'
dlx4_2='MKLSVLPPRSLLAPYTVLCCPPDSEKPRLSPELSERRPQAPAKKLRKPRTIYSSLQLQHLNQRFQHTQYLALPERAQLAAQLGLTQTQVKIWFQNKRSKYKKLLKQNSGGQEGDFPGRTFSVSPCSPPLPSLWDLPKAGTLPTSGYGNSFGAWYQHHSSDVLASPQMM'
pknox1_1='MMATQTLSIDSYQDGQQMQVVTELKTEQDPNCSEPDAEGVSPPPVESQTPMDVDKQAIYRHPLFPLLALLFEKCEQSTQGSEGTTSASFDVDIENFVRKQEKEGKPFFCEDPETDNLMVKAIQVLRIHLLELEKVNELCKDFCSRYIACLKTKMNSETLLSGEPGSPYSPVQSQQIQSAITGTISPQGIVVPASALQQGNVAMATVAGGTVYQPVTVVTPQGQVVTQTLSPGTIRIQNSQLQLQLNQDLSILHQDDGSSKNKRGVLPKHATNVMRSWLFQHIGHPYPTEDEKKQIAAQTNLTLLQVNNWFINARRRILQPMLDSSCSETPKTKKKTAQNRPVQRFWPDSIASGVAQPPPSELTMSEGAVVTITTPVNMNVDSLQSLSSDGATLAVQQVMMAGQSEDESVDSTEEDAGALAPAHISGLVLENSDSLQ'
pknox1_3='MVKAIQVLRIHLLELEKVNELCKDFCSRYIACLKTKMNSETLLSGEPGSPYSPVQSQQIQSAITGTISPQGIVVPASALQQGNVAMATVAGGTVYQPVTVVTPQGQVVTQTLSPGTIRIQNSQLQLQLNQDLSILHQDDGSSKNKRGVLPKHATNVMRSWLFQHIGHPYPTEDEKKQIAAQTNLTLLQVNNWFINARRRILQPMLDSSCSETPKTKKKTAQNRPVQRFWPDSIASGVAQPPPSELTMSEGAVVTITTPVNMNVDSLQSLSSDGATLAVQQVMMAGQSEDESVDSTEEDAGALAPAHISGLVLENSDSLQ'

# make regions dict for each. 
# As defined by disorder prediction or disorder prediction + small changes by Stephen. 
dlx4_1_regions={'idr_1':[0, 123], 'folded_1':[123, 180], 'idr_2':[180, len(dlx4_1)-1]}
dlx4_2_regions={'idr_1':[0, 51], 'folded_1':[51, 108], 'idr_2':[108, len(dlx4_2)-1]}
pknox1_1_regions={'idr_1':[0, 108], 'folded_1':[108, 155], 'idr_2':[155, 265], 'folded_3':[265, 325], 'idr_3':[325, len(pknox1_1)-1]}
pknox1_3_regions={'folded_1':[0,38], 'idr_1':[38, 148], 'folded_3':[148,208], 'idr_2':[208, len(pknox1_3)-1]}

# set loc of folder to save dodo files
dodo_loc='../dodo_structures'

# make dodo structures
pdb_from_pdb('dlx4_1.pdb', out_path=f'{dodo_loc}/dlx4_1_dodo.pdb', regions_dict=dlx4_1_regions, linear_placement=True)
pdb_from_pdb('dlx4_2.pdb', out_path=f'{dodo_loc}/dlx4_2_dodo.pdb', regions_dict=dlx4_2_regions, linear_placement=True)
pdb_from_pdb('pknox1_1.pdb', out_path=f'{dodo_loc}/pknox1_1_dodo.pdb', regions_dict=pknox1_1_regions, linear_placement=True)
pdb_from_pdb('pknox1_3.pdb', out_path=f'{dodo_loc}/pknox1_3_dodo.pdb', regions_dict=pknox1_3_regions, linear_placement=True)

