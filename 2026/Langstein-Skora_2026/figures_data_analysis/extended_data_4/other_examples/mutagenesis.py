# About
# --------------------------------------------------------------------------
#
# This module implements a mutagenesis protocol to generate new protein sequences
# that maintain similar interaction profiles to a given input sequence when 
# interacting with a fixed partner sequence. The protocol uses either the Mpipi
# or CALVADOS coarse-grained models to evaluate the interaction matrices of the
# sequences. The mutagenesis process involves generating libraries of mutants
# through random mutations and selecting those that closely resemble the wild type
# interaction profile based on a defined threshold. The module allows for multiple
# rounds of mutagenesis to accumulate mutations while preserving interaction characteristics.
#
# We note we are preparing a manuscript describing this protocol in detail (Flynn et al. in prep);
# if you feel this type of analysis is useful for your work, please reach out to us for our updated
# version of this code!
#
# ~alex (Nov 2025)


from finches import CALVADOS_frontend, Mpipi_frontend
import numpy as np

# nb - housetools is a holehouselab specific package; IF you want to run this code yourself we 
# recommend you contact us so we can provide you with a newer version of this code!
from housetools import nucleic_acid_tools 

# --------------------------------------------------------------------------
#
#
def compare(s1,s2, verbose=True):
    """
    Compare two sequences and count the number of differences between them.

    Parameters:
    -----------
    s1: str
        The first sequence to compare.

    s2: str
        The second sequence to compare.

    verbose: bool
        Whether to print the differences between the sequences.

    Returns:
    --------
    int
        The number of differences between the two sequences.
    """
    count_diffs = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count_diffs = count_diffs + 1
            if verbose:
                print(f"{i+1}: {s1[i]} vs. {s2[i]}")

    if verbose:
        print(f"{count_diffs} changes")
    return count_diffs

# --------------------------------------------------------------------------
#
#
def run_mutagenesis(input_sequence,
                    partner_sequence,
                    num_mutants,
                    delta_threshold=0.5,
                    num_loops=20,
                    target_size=30,
                    number_of_events=3,
                    mode='mpipi',
                    selection_type='maxdiff',
                    verbose=True):
    """

    Parameters:
    -----------

    input_sequence : str
        The input sequence that we want to mutate

    partner_sequence : str
        The partner sequence that encodes the interaction constraint
        on mutation. This sequence remains fixed during mutagenesis.

    num_mutants : int
        The number of mutants we want to generate.

    delta_threshold : float
        The threshold for the difference between the interaction matrix of
        the wild type and the mutant. If the difference is less than this
        threshold, we consider the mutant to be sufficiently close to the
        wild type.

    num_loops : int
        The number of loops we want to run the mutagenesis for. In each loop,
        we generate a library of mutants, and select one that is closest
        to the wild type. Number of loops is number of rounds of sequential
        mutagenesis.

    target_size : int
        The size of the library of mutants we want to generate in each round.
        If we're trying to obtain mutants with many mutations we might want
        to increase this number to increase the chances of finding a mutant
        that is close to the wild type.

    number_of_events : int
        The number of mutations we want to make in each round of mutagenesis.
        When larger this becomes harder to find a mutant that is close to the
        wild type, but we can generate mutants with more mutations.

    mode : str
        The mode of the model, either 'mpipi' or 'calvados'.

    selection_type : str
        The type of selection we want to use to. Options are 
       
        - 'maxdiff' : delta reflects the max diff between any pair of 
                        elements in the interaction matrix

        - 'sumdiff' : delta reflects the sum of the differences between
                        the interaction matrix of the wild type and the
                        mutant

        - 'attractive_sumdiff' : delta reflects the sum of the differences
                        between the interaction matrix of the wild type and
                        the mutant, but only for the attractive interactions.


    verbose: bool
        Whether to print information during mutagenesis.

    Returns:
    --------
    dict
        A dictionary of the new sequences that were generated during mutagenesis.

    """


    # define model
    if mode == 'mpipi':
        xf = Mpipi_frontend()
    elif mode == 'calvados':
        xf = CALVADOS_frontend()
    else:
        raise ValueError('mode must be either mpipi or calvados')
    
    # get the interaction matrix for the input sequence
    wt = xf.intermolecular_idr_matrix(input_sequence, partner_sequence)

    # initialize the new sequence library
    new_seq_library = {}

    # cycle through the number of mutants
    for idx, _ in enumerate(range(num_mutants)):

        # initialize the start sequence
        start_seq = input_sequence

        # cycle through the number of loops
        for l in range(num_loops):

            # initialize an empty start sequence
            new_start_seq = ''

            # build a mutant library target size = number of mutants we try and generate, but often we don't get that many
            # number of events = number of mutations we try and make
            new_seqs = nucleic_acid_tools.build_mutant_library(start_seq, target_size=target_size, number_of_events=number_of_events)
                                                               
            # cycle through the new sequences
            for s in new_seqs:

                # get the interaction matrix for the new sequence
                x = xf.intermolecular_idr_matrix(s, partner_sequence)

                if selection_type == 'maxdiff':
                    # calculate the difference between the interaction matrix of the new sequence and the wild type, where
                    # we ask which is the pairwise difference that's largest
                    delta = np.max(abs(x[0][0]- wt[0][0]))

                elif selection_type == 'sumdiff':

                    delta = abs(np.sum(x[0][0]) - np.sum(wt[0][0]))

                elif selection_type == 'attractive_sumdiff':
                    delta = abs(np.sum(x[0][0][x[0][0] < 0]) - np.sum(wt[0][0][wt[0][0] < 0]))

                else:
                    raise ValueError('selection_type invalid')
                
                # if we find a sequence that's sufficiently close to the wild type, we break
                if delta < delta_threshold:
                    new_start_seq = s
                    break

            # new_start_seq was initialized as an empty string, so if we didn't find a sequence that was close enough
            # to the wild type, we don't update the start sequence, otherwise we do, so we perform sequential rounds of
            # mutagenesis on the same sequence
            if len(new_start_seq) > 0:
                start_seq = new_start_seq

        # FINALLY, check to see if we found a sequence that was close enough to the wild type in the final
        # loop, if not, we just use the sequence from the last loop
        if len(new_start_seq) > 0:
            new_seq = new_start_seq
        else:
            new_seq = start_seq

        # count the number of mutations
        d = compare(input_sequence, new_seq, verbose=False)

        # add the new sequence to the library
        new_seq_library[f"s_{idx}"] = new_seq

        # determine the average number of mutations
        all_comps = []
        for k in new_seq_library:
            all_comps.append(compare(input_sequence, new_seq_library[k], verbose=False))

        if verbose:
            print(f'Sequence {idx} has {d} mutations (running average: {round(np.mean(all_comps),2)})')

    return new_seq_library
            

