import protfasta
import metapredict as meta



########################################################################################
#
class ProteinException(Exception):
    """
    Simple exception class 
    """
    pass
########################################################################################


########################################################################################
#
class Protein:
    """
    Simple class that functions as a datastructure to hold protein sequence and 
    disorder information
    """

    def __init__(self, sequence, disorder_binary, disorder_score=None):


        # sanity check
        if len(sequence) != len(disorder_binary):
            raise  ProteinException('Sequence and disorder binary scores did not match')

        # sanity check
        if disorder_score is not None:
            if len(sequence) != len(disorder_score):
                raise  ProteinException('Sequence and disorder score scores did not match')

        self.__sequence = sequence
        self.__disorder_score = disorder_score
        self.__disorder_binary = disorder_binary        

    @property
    def sequence(self):
        return self.__sequence

    @property
    def disorder_score(self):
        if self.__disorder_score is None:
            return [-1.0]*len(self.__sequence)
        else:
            return self.__disorder_score

    @property
    def disorder_binary(self):
        return self.__disorder_binary


    def __str__(self):
        return f"Protein with {len(self.sequence)} residues"

    def __repr__(self):
        return str(self)

########################################################################################



# ..................................................................................
#
#
def parse_file(fn):
    """
    Function that parses a CAID disorder prediction file and returns a dictionary of key-value pairs
    where keys are DISPROT IDs and values are Protein objects .

    Expectation is that CAID disorder prediction files have the format
    
    >DISPROT_ID_A
    <res_id1> <residue1> <disorder_score1> <disorder_binary1>
    <res_id2> <residue2> <disorder_score2> <disorder_binary2>
    <res_id3> <residue3> <disorder_score3> <disorder_binary3>

    ...

    >DISPROT_ID_B
    <res_id1> <residue1> <disorder_score1> <disorder_binary1>
    <res_id2> <residue2> <disorder_score2> <disorder_binary2>
    <res_id3> <residue3> <disorder_score3> <disorder_binary3>

    ...

    Although there need not necessarily be four columns - so predictors offer only a binary score
    as opposed to a binary and a continous score.


    Parameters
    -------------

    fn : str
        Filename for a CIAD disorder rediction file 

    """

    # initialize the protein dictionary, which will ultimately be the dictionary we return
    protein_dict = {}

    # read file in
    with open(fn, 'r') as fh:
        content = fh.readlines()
    
    # grab first disprot ID from the data
    current_pid = content[0][1:].strip()
    current_seq = ''
    current_score = []
    current_score_bin = []

    # cycle through each line in the file starting from the 2nd line
    for line in content[1:]:
        
        ## The if statement below happens when we find a new protein entry (i.e. a new >DID line)
        # if we find an ID line
        if line[0] == '>':
            
            # if we'd previously found continous scores then use the scores
            if len(current_score) > 0: 
                protein_dict[current_pid] = Protein(current_seq, current_score_bin, current_score)
            else:
                protein_dict[current_pid] = Protein(current_seq, current_score_bin)
            
            # update and reset the various values
            current_pid = line[1:].strip()
            current_seq = ''
            current_score = []
            current_score_bin = []
            
        else:
            # split the line and extract out the relevant params
            sl = line.strip().split()
            
            if len(sl) == 4:
                current_seq = current_seq + sl[1]
                current_score.append(float(sl[2]))
                current_score_bin.append(int(sl[3]))

            elif len(sl) == 3:
                current_seq = current_seq + sl[1]
                current_score_bin.append(int(sl[2]))
    
    # finally, we have a final entry in the file which we have to parse (same as done in the line[0] == ">": 
    # if statement
    if len(current_score) > 0: 
        protein_dict[current_pid] = Protein(current_seq, current_score_bin, current_score)
    else:
        protein_dict[current_pid] = Protein(current_seq, current_score_bin)

    return protein_dict



# ..................................................................................
#
#
def parse_disprot(verbose=False, datadir='data'):
    """
    Function that parses the disprot database. Note this reads three files DISPROT_* filenames defined as
    constant variables in this file at the top. As such, it doesn't actually take any arguments but is 
    obviously not super transferrable as written...

    Importantly, this function divides each residue associated with a  DISPROT sequence 
    into one of three possible annotations:

    0  - structured, unambigiously identified as structured from PDB data and not disordered
    1  - disordered, unamvigiously identified as disordered and NOT structured based on PDB data
    -2 - ambigious, residue for which EITHER there's conflicting PDB and disorder data OR no data from either

    Any DISPROT entries that lack one (or more) of sequence, disorder, and structural data are ignored and a 
    warning is printed.

    Returns
    ---------
    dict
        Returns a dictionary of key-value pairs where keys are DISPROT IDs and values are Protein objects, 
        which possess sequence and disorder data defined as described above.


    """

    DISPROT_SEQ = f'{datadir}/disprot-2018-11-seq.fasta'
    DISPROT_DIS = f'{datadir}/disprot-2018-11-disorder.fasta'
    DISPROT_PDB = f'{datadir}/pdb-atleast.fasta'


    # ----------------------------------------
    def convert_to_list(ls):
        """
        Internal function that converts a string to a list
        where string entries that are '-' are set to 0 while
        anything else is set to 1

        Parameters
        ------------
        ls : str
            Sequence string 

        Returns
        ------------
        list
            Returns a binary list of 0s and 1s


        """
        x = []
        for k in ls:
            if k == '-':
                x.append(0)
            else:
                x.append(1)
        return x
    # ----------------------------------------

    # define the Disprot database
    database = {}

    ## Extract out disprot sequence data (i.e. amino acid sequences. Note this converts non-standard
    # amino acids to standard AAs (of which there's only 1-2 examples across the dataset)
    def hp(s):
        return(s.strip().split()[0])
    dis_seq = protfasta.read_fasta(DISPROT_SEQ, header_parser=hp, invalid_sequence_action='convert')


    ## Extract out disprot disorder data (where residues are either '-' ("not disordered") or 'D' ("disordered")
    def hp(s):
    
        if s == 'this test string should work':
            return s
        did = s.split('|')[1]                    
        return(did)
    dis_disorder = protfasta.read_fasta(DISPROT_DIS, header_parser=hp, invalid_sequence_action='ignore')    

    ## Extract out disprot structural data (where residues are either '-' ("not structured") or 'S' ("structured")
    def hp(s):
        return(s.strip().split()[0])
    dis_pdb = protfasta.read_fasta(DISPROT_PDB, header_parser=hp, invalid_sequence_action='ignore')    

    # i here is a DisProt ID and a key into each of the sequence dictionaries
    for i in dis_disorder:

        # if we found an ID that's not found in the disprot disorder or structural data then we just
        # skip over it
        if i not in dis_seq: 
            if verbose:
                print(f'Warning - could not find sequence for {i}')
        elif i not in dis_pdb:
            if verbose:
                print(f'Warning - could not find structural data for for {i}')
        else:

            # get the sequence
            seq = dis_seq[i]

            # define the initial dis list, which is a 0 or 1 based on if a residue is disorder (1) or 
            # not
            dis = convert_to_list(dis_disorder[i])

            ## this for-loop updates the dis list such that each position is set to 
            #
            #  0 => structured (unambigious)
            #  1 => disordered (unambigious)
            # -2 => Amigious 

            # Specifically, we cycle over each residue position
            for residue_position in range(len(dis)):

                # if this residue is structured (based on PDB data)
                if dis_pdb[i][residue_position] == 'S':

                    # AND if the residue was NOT disordered then we take this as unambigiously structured
                    if dis[residue_position] == 0:
                        pass
                        
                    # but if the disorder data said this was disordered we assign it to -2 (unsure)
                    else:
                        dis[residue_position] = -2

                # alternatively, if this residue is NOT structured
                else:

                    # and if this residue is disordered, then we take this as unamvigiously disordered
                    if dis[residue_position] == 1:
                        pass

                    # but if it wasn't assigned to be disordered we say it's unsure
                    else:
                        dis[residue_position] = -2
                        
                        
                    
            # assign this entry as  protein object to the database
            database[i] = Protein(seq, dis)


    if verbose:
        print(f'Parsing complete - DISPROT dataset contains {len(database)} protein entries')

    return database

    


# ..................................................................................
#
#
def extract_binary_disorder(seq, 
                            threshold=None, 
                            mode='legacy'):
    """
    This function takes an amino acid sequence and returns a binary disorder list
    where each residue is assigned either disordered or ordered as defined by
    metapredicts's predict_disorder_domains() function.

    Parameters
    ------------
    seq : str
        Amino acid sequence

    threshold : float
        Threshold value passed to metapredicts disorder prediction algorithm. 
        If not provided default values for the methods are used

    mode : str
        Selector for the method used. Valid optoins are:

            legacy (original metapredict)

            v2 (new hybrid method that combines metapredict and ppLDDT)

        Default = 'legacy'

    Parameters
    -------------
    list
        List of 0s and 1s where disorder =1 and not disordered = 0.
    

    """
    
    if mode == 'legacy':
        f = meta.predict_disorder_domains(seq, disorder_threshold=threshold, legacy=True)
    elif mode == 'v2':        
        f = meta.predict_disorder_domains(seq, disorder_threshold=threshold)
    else:
        raise Exception('Modes must be legacy or v2')

    # get IDR boundaries
    idrs = f.disordered_domain_boundaries
    binary_scores = []

    # for each position in the sequence
    for i in range(len(seq)):

        # set found to False
        found = False       

        # for each set of disorder boundaries identified by metapredict ask if
        # the curret position is within the boundaries- if yes set found
        # to True
        for d in idrs:
            if i >= d[0] and i < d[1]:
                binary_scores.append(1)
                found =True
                break
                
        # if we did not find a predicted IDR that position $i falls into 
        # then set this to 9
        if found is False:
            binary_scores.append(0)

    return binary_scores
        

# ..................................................................................
#
#
def calculate_residue_count(truth, prediction):
    """
    Calculate count of residues in a prediction ($prediction) that are right or wrong based on disprot data 
    $(truth).

    Note we IGNORE residues where disport is ambigious, so 

    Parameters
    ------------
    truth : list
        List where each position is either 0 (structured), 1 (disordered), or -2 (ambigous)
    
    prediction : list
        List where ach position is either predicted to be structured (0) or disordered (1)

    Returns
    -----------
    (valid, invalid) tuple
        Returns a tuple where the FIRST entry is number of residues where prediction and truth are
        the same, and invalid is the number of residues where prediction and truth were different, 
        ignoring any positions where ground truth was ambigious.
        
    """

    valid = 0
    invalid = 0

    for i in range(len(truth)):

        # we take advanatge of the fact that the prediction
        # can only be 0 or 1, so if they're the same we take
        # this to be a true positive prediction
        if truth[i] == prediction[i]:
            valid = valid + 1

        # if truth is ambigious then we don't count this as 
        # anything
        elif truth[i] == -2:
            pass

        # HOWEVER, if we get here the truth was unamvigiously *something* and 
        # that something was not what the prediction expect, which we take as a
        # a failure
        else:
            invalid = invalid + 1

    return (valid, invalid)
            
        
