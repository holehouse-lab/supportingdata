"""Per-residue accuracy (correctly scored residues per 100) on CAID3 Disorder-PDB.

The reader functions are taken from the original Figure S7 notebook (which scored
CAID2). The changes are: the ground truth is passed explicitly, the CAID3 scoring
rules match Figure 1C (each method's own binary calls, or a fixed cutoff for
raw-score-only methods), and metapredict is scored from metapredict 3.1.0
predictions rather than the files distributed by CAID.
"""

from __future__ import annotations

import csv

from figtools import metapredict_caid_path, shared_data_path

# Raw-score-only CAID3 methods are thresholded at 0.5 (every CAID3 raw score is on
# a [0, 1] scale); the AlphaFold RSA / pLDDT tracks reuse the CAID2 thresholds.
CAID3_CUTOFF_OVERRIDES = {
    "AlphaFold-rsa": 0.581,
    "AlphaFold3-rsa": 0.581,
    "AlphaFold-pLDDT": 0.312,
    "AlphaFold3-pLDDT": 0.312,
}
DEFAULT_CUTOFF = 0.5

METAPREDICT_NAMES = {
    "Metapredict-v1": "metapredict V1",
    "Metapredict-v2": "metapredict V2",
    "Metapredict-v3": "metapredict V3",
}


def read_caid_seq_disorder(caid_fasta_file: str) -> dict:
    """Read a CAID reference (header, sequence, label string per entry).

    Parameters
    ----------
    caid_fasta_file : str
        Path to the reference; labels are ``1`` disordered, ``0`` ordered, ``-``
        not assessed.

    Returns
    -------
    dict
        ``{'>ID': {'sequence': str, 'scores': str}}``.
    """
    caid_seqs_scores = {}
    with open(caid_fasta_file, "r") as fh:
        lines = fh.readlines()
    for n in range(0, len(lines), 3):
        caid_seqs_scores[lines[n].strip()] = {"sequence": lines[n + 1].strip(),
                                              "scores": lines[n + 2].strip()}
    return caid_seqs_scores


def read_caid_convert_get_predictor_vals(caid_file: str, ref_data: dict) -> list:
    """Binary predictions and labels for the assessed residues, from submitted calls.

    Parameters
    ----------
    caid_file : str
        ``.caid`` file with a populated 4th (binary) column.
    ref_data : dict
        Ground truth from :func:`read_caid_seq_disorder`.

    Returns
    -------
    list
        ``[predictor_vals, caid_vals]``, matched lists of 0/1 values.
    """
    caid_vals = []
    predictor_vals = []
    with open(caid_file, "r") as fh:
        preds = fh.read().split(">")
    for p in preds:
        if p == "":
            continue
        pred_lines = p.split("\n")
        name = pred_lines[0].split("\t")[0]
        curdat = ref_data[f">{name}"]["scores"]
        score_ind = 0
        for v in pred_lines:
            splitvals = v.split("\t")
            if len(splitvals) == 4:
                actual = curdat[score_ind]
                if actual != "-":
                    caid_vals.append(int(actual))
                    predictor_vals.append(int(splitvals[-1]))
                score_ind += 1
    return [predictor_vals, caid_vals]


def read_caid_convert_get_vals_allow_altered_cutoff(caid_file: str, ref_data: dict,
                                                    cutoff: float = 0.5) -> list:
    """Binary predictions and labels, thresholding raw scores where no call exists.

    Parameters
    ----------
    caid_file : str
        ``.caid`` file (3 or 4 columns; the binary column may be blank).
    ref_data : dict
        Ground truth from :func:`read_caid_seq_disorder`.
    cutoff : float, optional
        Disorder threshold applied to raw scores (default 0.5).

    Returns
    -------
    list
        ``[predictor_vals, caid_vals]``, matched lists of 0/1 values.
    """
    caid_vals = []
    predictor_vals = []
    with open(caid_file, "r") as fh:
        preds = fh.read().split(">")
    for p in preds:
        if p == "":
            continue
        pred_lines = p.split("\n")
        name = pred_lines[0].split("\t")[0]
        curdat = ref_data[f">{name}"]["scores"]
        score_ind = 0
        for v in pred_lines:
            splitvals = v.split("\t")
            if len(splitvals) < 3:
                continue
            if len(splitvals) == 4 and splitvals[-1] != "":
                curpred = int(splitvals[-1])
            else:
                rawval = float(splitvals[2])
                curpred = 1 if rawval >= cutoff else 0
            actual = curdat[score_ind]
            if actual != "-":
                caid_vals.append(int(actual))
                predictor_vals.append(curpred)
            score_ind += 1
    return [predictor_vals, caid_vals]


def caid3_has_binarized_scores(caid_file: str) -> bool:
    """True if a ``.caid`` file provides a usable binary column (4th field populated)."""
    with open(caid_file, "r") as fh:
        for line in fh:
            splitvals = line.rstrip("\n").split("\t")
            if len(splitvals) == 4 and splitvals[-1] != "":
                return True
    return False


def correct_per_100(predictor_vals: list, caid_vals: list) -> float:
    """Number of correctly labelled residues per 100 assessed residues."""
    assert len(predictor_vals) == len(caid_vals), "prediction/label length mismatch"
    correct = sum(1 for pred, actual in zip(predictor_vals, caid_vals) if pred == actual)
    return 100 * correct / len(caid_vals)


def score_caid3_field(reference_fasta: str, results_csv: str) -> dict:
    """Correct residues per 100 for every CAID3 disorder predictor.

    Parameters
    ----------
    reference_fasta : str
        CAID3 Disorder-PDB reference.
    results_csv : str
        CAID3 disorder results table; only its ``predictor`` column (the list of
        disorder predictors, binding predictors excluded) is used.

    Returns
    -------
    dict
        Predictor label -> correct residues per 100 (metapredict from 3.1.0).
    """
    ref_data = read_caid_seq_disorder(reference_fasta)
    with open(results_csv) as fh:
        names = [row["predictor"] for row in csv.DictReader(fh)]

    scores = {}
    for name in names:
        if name in METAPREDICT_NAMES:
            caid_file = metapredict_caid_path("caid3", name[-2:].upper())
        else:
            caid_file = shared_data_path("caid3", f"{name}.caid")
        if caid3_has_binarized_scores(caid_file):
            vals = read_caid_convert_get_predictor_vals(caid_file, ref_data)
        else:
            cutoff = CAID3_CUTOFF_OVERRIDES.get(name, DEFAULT_CUTOFF)
            vals = read_caid_convert_get_vals_allow_altered_cutoff(caid_file, ref_data, cutoff)
        scores[METAPREDICT_NAMES.get(name, name)] = correct_per_100(*vals)
    return dict(sorted(scores.items(), key=lambda item: item[1], reverse=True))
