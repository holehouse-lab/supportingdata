"""Code used to build the 'combined disorder scores' that metapredict V3 was trained on.

``moving_average`` and ``meta_predict_hybrid_v3`` are copied verbatim from
``metapredict/analysis/creating_V2_and_V3_scores.py`` in the metapredict repository
(https://github.com/idptools/metapredict). That module is kept in the repository
for reproducibility but is not shipped in the installed package (metapredict
>= 3.1.0), so the two functions are reproduced here unchanged. The module's
``alphaPredict`` import is omitted because neither function uses it.
"""

# import local modules
from metapredict.backend.predictor import predict

# import everything else
import numpy as np
from scipy.signal import savgol_filter
from tqdm import tqdm


def moving_average(scores, window_size):
    """
    Smooths an array of scores over a user-specified sliding window size.

    Parameters
    ----------
    scores : array
        Array of scores to be smoothed.
    window_size : int
        Size of the sliding window for smoothing.

    Returns
    -------
    smoothed_scores : ndarray
        Smoothed scores with the same length as the input array.
    """
    if window_size % 2 == 0:
        raise ValueError("Window size must be an odd number.")

    half_window = window_size // 2
    pad_width = ((half_window, half_window),)  # Pad equally on both sides
    padded_scores = np.pad(scores, pad_width, mode='reflect')
    smoothed_scores = np.convolve(padded_scores, np.ones(window_size) / window_size, mode='valid')
    return smoothed_scores



def meta_predict_hybrid_v3(inputs, metapredict_version='v1', vmax_cutoff=0.5,
    base=0.35, top=0.95, windowed_averaging=25):
    """
    Code used to generate the scores used to train metapredict V3.
    All default values in the function are as used to make V3 scores.
    Real pLDDT scores (not predicted) were used to make V3.

    Parameters
    ------------
    inputs : dict
        Sequences as the keys and the plddt scores as the values.

    metapredict_version : str
        The version of metapredict to use.
        Default = 'v1'

    vmax_cutoff : float
        The cutoff value for the maximum value.
        Default = 0.5

    base : float
        The base value for the plddt scores.
        Default = 0.35

    top : float
        The top value for the plddt scores.
        Default = 0.95

    windowed_averaging : Int
        Whether to smooth scores. Default is 25.
        Value to make V3 scores for training was 25.

    Returns
    ----------
    dict of nd.nparrays
        Returns a dict of arrays with the sequence as the key and the
        scores as the values.
    """
    # make dict to hold sequences and plddt / disorder scores
    seq_plddt_disorder_dict={}

    # make sure a dict is input
    if isinstance(inputs, dict):
        sequences=list(inputs.keys())
        disorder_scores = predict(sequences, version=metapredict_version, round_values=False)
        for n in range(0, len(disorder_scores)):
            seq = disorder_scores[n][0]
            disorder = disorder_scores[n][1]
            plddt=np.array(inputs[seq])
            seq_plddt_disorder_dict[seq] = {'disorder':disorder,'pLDDT': plddt}
    else:
        raise Exception('only a dict is accepted as input for this function.')

    # now we need to iterate over the dictionary and get hybrid scores.
    hybrid_scores = {}
    for seq in tqdm(seq_plddt_disorder_dict):
        pLDDT = seq_plddt_disorder_dict[seq]['pLDDT']
        disorder = seq_plddt_disorder_dict[seq]['disorder']

        # defines functionally the limits that the predicted pLDDT score
        # can be between
        base = base
        top  = top

        # this normalizes so the pLDDT score ends up being renormalized
        # to be between 0 and 1, converted into an effective disorder score
        d = (pLDDT - base)*(1/(top-base))

        # means value of 1 = disordered and 0 is disordered
        d = 1 - d

        hybrid_vals = []
        for idx in range(len(d)):

            vmax = max([disorder[idx], d[idx]])

            # if the largest disorder score for either is above 0.5 go
            # with that score as the consensus
            if vmax > vmax_cutoff:
                hybrid_vals.append(vmax)

            # else go with the smallest disorder score as the consensus
            else:
                hybrid_vals.append(min([disorder[idx], d[idx]]))

        # smooth with a window size of 7 to kill discontinuities
        smoothed = savgol_filter(hybrid_vals,7,3)

        # flatten so all values in the bounds of 0 to 1
        smoothed = np.where(smoothed<0, 0, smoothed)
        smoothed = np.where(smoothed>1, 1, smoothed)

        # moving averaging
        if windowed_averaging != None:
            smoothed=moving_average(smoothed, windowed_averaging)

        hybrid_scores[seq]=np.round(smoothed,4).tolist()
    return hybrid_scores
