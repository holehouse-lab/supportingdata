"""CAID reading, MCC and ranking-panel helpers for Figure 1A-C.

These functions are taken from the original ``Figure_1A_1B.ipynb`` notebook. The
changes are: (i) the reference reader opens the path it is given (the original
re-joined the path onto its own directory, which breaks for files in ``data/``),
(ii) the ground-truth dictionary is passed explicitly rather than defaulting to a
CAID2 global, (iii) the ranking panel follows the manuscript style (solid black
bars, dataset name as an in-axis label rather than a title, shared y-limits), and
(iv) metapredict is always scored from metapredict 3.1.0 predictions
(:func:`metapredict_mcc`) rather than hard-coded values or the CAID-distributed
files.
"""

from __future__ import annotations

import math

import numpy as np

from figtools import FIELD_COLOR, V1_COLOR, V2_COLOR, V3_COLOR, metapredict_caid_path

# Panels 1A-C share one MCC axis so the three CAID rounds can be compared by eye.
# Upper limit chosen to clear the best predictor in any round (max MCC ~0.80);
# plot_mcc_panel asserts that nothing is clipped.
SHARED_MCC_YLIM = (0.0, 0.84)

# CAID3 methods that submitted only raw scores need a disorder cutoff. Every CAID3
# raw score is on a [0, 1] scale, so 0.5 is the default; the AlphaFold RSA / pLDDT
# tracks are on physical scales and reuse the CAID2 thresholds (CAID2, Conte et al.
# 2023).
CAID3_CUTOFF_OVERRIDES = {
    "AlphaFold-rsa": 0.581,
    "AlphaFold3-rsa": 0.581,
    "AlphaFold-pLDDT": 0.312,
    "AlphaFold3-pLDDT": 0.312,
}
DEFAULT_CUTOFF = 0.5

# CAID3 file names for the metapredict versions -> the labels used in every panel.
CAID3_METAPREDICT_NAMES = {
    "Metapredict-v1": "metapredict V1",
    "Metapredict-v2": "metapredict V2",
    "Metapredict-v3": "metapredict V3",
}


def read_caid_seq_disorder(caid_file: str) -> dict:
    """Read a CAID reference FASTA (header, sequence, label string per entry).

    Parameters
    ----------
    caid_file : str
        Path to the reference file. Labels are ``1`` disordered, ``0`` ordered and
        ``-`` not assessed.

    Returns
    -------
    dict
        ``{'>ID': {'sequence': str, 'scores': str}}``.
    """
    caid_seqs_scores = {}
    with open(caid_file, "r") as fh:
        lines = fh.readlines()
    for n in range(0, len(lines), 3):
        caid_seqs_scores[lines[n].strip()] = {
            "sequence": lines[n + 1].strip(),
            "scores": lines[n + 2].strip(),
        }
    return caid_seqs_scores


def calc_mcc(predictor_scores: list, actual_scores: list) -> float:
    """Matthews correlation coefficient from binary predictions and labels.

    Parameters
    ----------
    predictor_scores : list of int
        Predicted labels, disorder = 1, order = 0.
    actual_scores : list of int
        Ground-truth labels, disorder = 1, order = 0.

    Returns
    -------
    float
        MCC as defined in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4406047/
    """
    true_positive = 0
    false_positive = 0
    true_negative = 0
    false_negative = 0
    for i in range(0, len(predictor_scores)):
        predicted = predictor_scores[i]
        actual = actual_scores[i]
        if predicted == 1 and actual == 1:
            true_positive = true_positive + 1
        elif predicted == 1 and actual == 0:
            false_positive = false_positive + 1
        elif predicted == 0 and actual == 0:
            true_negative = true_negative + 1
        elif predicted == 0 and actual == 1:
            false_negative = false_negative + 1

    numerator = (true_positive * true_negative) - (false_positive * false_negative)
    denomenator_before_square = (
        (true_positive + false_positive)
        * (true_negative + false_positive)
        * (true_positive + false_negative)
        * (true_negative + false_negative)
    )
    denomenator = math.sqrt(denomenator_before_square)
    return numerator / denomenator


def read_caid_convert_get_MCC_allow_altered_cutoff(caid_file: str, ref_data: dict,
                                                    cutoff: float = 0.5) -> float:
    """MCC for a ``.caid`` file, thresholding raw scores where no binary call exists.

    CAID files are inconsistently formatted between methods: lines are
    ``index, residue, raw_score[, binary_call]`` and either field can be blank. A
    submitted binary call is used when present; otherwise the raw score is
    thresholded at ``cutoff``.

    Parameters
    ----------
    caid_file : str
        Path to the ``.caid`` prediction file.
    ref_data : dict
        Ground truth from :func:`read_caid_seq_disorder`.
    cutoff : float, optional
        Disorder threshold applied to raw scores (default 0.5).

    Returns
    -------
    float
        Dataset-wide (pooled per-residue) MCC.
    """
    caid_vals = []
    predictor_vals = []
    with open(caid_file, "r") as fh:
        preds = fh.read().split(">")
    for p in preds:
        if p == "":
            continue
        pred_lines = p.split("\n")
        name = pred_lines[0]
        if "\t" in name:
            name = name.split("\t")[0]
        curdat = ref_data[f">{name}"]["scores"]
        score_ind = 0
        for v in pred_lines:
            splitvals = v.split("\t")
            if len(splitvals) < 3:
                continue
            if len(splitvals) == 4:
                if splitvals[-1] == "":
                    curpred = 1 if float(splitvals[-2]) >= cutoff else 0
                else:
                    curpred = int(splitvals[-1])
            else:
                curpred = 1 if float(splitvals[-1]) >= cutoff else 0
            actual = curdat[score_ind]
            if actual != "-":
                caid_vals.append(int(actual))
                predictor_vals.append(curpred)
            score_ind += 1
    return calc_mcc(predictor_vals, caid_vals)


def read_caid_convert_get_MCC(caid_file: str, ref_data: dict) -> float:
    """MCC for a ``.caid`` file using only the method's own submitted binary calls.

    Parameters
    ----------
    caid_file : str
        Path to the ``.caid`` prediction file (must have a 4th, binary column).
    ref_data : dict
        Ground truth from :func:`read_caid_seq_disorder`.

    Returns
    -------
    float
        Dataset-wide (pooled per-residue) MCC.
    """
    caid_vals = []
    predictor_vals = []
    with open(caid_file, "r") as fh:
        preds = fh.read().split(">")
    for p in preds:
        if p == "":
            continue
        pred_lines = p.split("\n")
        name = pred_lines[0]
        if "\t" in name:
            name = name.split("\t")[0]
        curdat = ref_data[f">{name}"]["scores"]
        score_ind = 0
        for v in pred_lines:
            splitvals = v.split("\t")
            if len(splitvals) == 4:
                curpred = int(splitvals[-1])
                actual = curdat[score_ind]
                if actual != "-":
                    caid_vals.append(int(actual))
                    predictor_vals.append(curpred)
                score_ind += 1
    return calc_mcc(predictor_vals, caid_vals)


def caid3_has_binarized_scores(caid_file: str) -> bool:
    """True if a ``.caid`` file provides a usable binary column (4th field populated)."""
    with open(caid_file, "r") as fh:
        for p in fh.read().split(">"):
            if p == "":
                continue
            for line in p.split("\n"):
                splitvals = line.split("\t")
                if len(splitvals) == 4 and splitvals[-1] != "":
                    return True
    return False


def metapredict_mcc(caid_round: str, ref_data: dict) -> dict:
    """MCC of metapredict V1, V2 and V3 (3.1.0 predictions) on one CAID round.

    Uses the binary calls in the 3.1.0 prediction files
    (``data/metapredict_caid/domains/<round>/``), where residues inside a predicted
    IDR are disordered.

    Parameters
    ----------
    caid_round : str
        ``"caid1"``, ``"caid2"`` or ``"caid3"``.
    ref_data : dict
        Ground truth for that round from :func:`read_caid_seq_disorder`.

    Returns
    -------
    dict
        ``{'metapredict V1': mcc, 'metapredict V2': mcc, 'metapredict V3': mcc}``.
    """
    return {f"metapredict {version}": read_caid_convert_get_MCC(
                metapredict_caid_path(caid_round, version), ref_data)
            for version in ["V1", "V2", "V3"]}


def plot_mcc_panel(ax, mcc_dict: dict, label: str, top_n: int | None = None,
                   ylim: tuple = SHARED_MCC_YLIM) -> None:
    """Draw a CAID-style MCC ranking onto ``ax``.

    Bars are ordered high to low (optionally truncated to ``top_n``), drawn in
    black, with any metapredict version recoloured in its signature colour. The
    dataset ``label`` and ``n`` (number of predictors shown, metapredict included)
    are written inside the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    mcc_dict : dict
        ``{predictor name: MCC}``; metapredict entries must be named
        ``'metapredict V1'`` / ``'V2'`` / ``'V3'``.
    label : str
        Dataset label written in the top-right corner.
    top_n : int or None, optional
        Keep only the ``top_n`` highest-MCC predictors.
    ylim : tuple of float, optional
        y-axis limits; defaults to the limits shared by panels 1A-C.
    """
    ordered = dict(sorted(mcc_dict.items(), key=lambda item: item[1], reverse=True))
    names = list(ordered.keys())
    scores = list(ordered.values())
    if top_n is not None:
        names = names[:top_n]
        scores = scores[:top_n]
    assert max(scores) < ylim[1], f"MCC {max(scores):.3f} would be clipped by ylim {ylim}"

    x = np.arange(1, len(scores) + 1)
    ax.bar(x, scores, color=FIELD_COLOR)
    for mp, mp_color in [("metapredict V3", V3_COLOR), ("metapredict V2", V2_COLOR),
                         ("metapredict V1", V1_COLOR)]:
        if mp in names:
            ax.bar(names.index(mp) + 1, scores[names.index(mp)], color=mp_color,
                   linewidth=0.5, edgecolor="black")

    ax.set_ylabel("MCC", fontsize=9)
    ax.set_xticks(x, names, fontsize=7, rotation=45, horizontalalignment="right",
                  rotation_mode="anchor")
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlim([0.5, len(scores) + 0.5])
    ax.set_ylim(ylim)
    ax.text(0.99, 0.96, f"{label}  ($n = {len(names)}$)", transform=ax.transAxes,
            ha="right", va="top", fontsize=8)
