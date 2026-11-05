"""Helpers for Figure S15: metapredict V3 on synthetic / designed proteins.

The plotting and fraction-disordered code is taken from the original Figure S8
notebooks (``hallucination_sequences``, ``pesce_sequences``, ``tretyachenko_analysis``
and ``goose_analysis``). Changes: every predicted IDR is shaded (the original loop
only kept the last IDR; identical output for these sequences, which have at most
one), traces are drawn in the metapredict V3 colour, and the protein label is
written inside the axes instead of being added at assembly.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import metapredict as meta
import numpy as np

from figtools import V3_COLOR

METAPREDICT_VERSION = 3


def plot_disorder(ax: plt.Axes, seq: str, label: str, start: int | None = None) -> None:
    """Plot a metapredict V3 disorder profile with predicted IDRs shaded.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    seq : str
        Protein sequence.
    label : str
        Identifier written just above the top-right corner of the axes.
    start : int or None, optional
        1-indexed residue of the native initiator methionine; if given, a dashed
        line marks where the N-terminal His-tag ends.
    """
    d = meta.predict_disorder_domains(seq, version=METAPREDICT_VERSION)
    idx = np.arange(1, len(seq) + 1)
    for idr_start, idr_end in d.disordered_domain_boundaries:
        idr_idx = np.arange(1 + idr_start, idr_end + 1)
        ax.fill_between(idr_idx, d.disorder[idr_start:idr_end], color=V3_COLOR, alpha=0.4, linewidth=0)
    ax.plot(idx, d.disorder, color=V3_COLOR, linewidth=0.8)
    if start is not None:
        ax.plot([start, start], [0, 1], '--k', linewidth=0.5)
    ax.set_xlim([idx[0], idx[-1]])
    ax.set_ylim([0, 1])
    ax.set_yticks([0, 1])
    ax.set_ylabel('Disorder', fontsize=8)
    ax.text(1.0, 1.02, label, transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    ax.tick_params(labelsize=7)


def fraction_disordered(seq: str) -> float:
    """Fraction of residues inside metapredict V3-defined IDRs.

    Parameters
    ----------
    seq : str
        Protein sequence.

    Returns
    -------
    float
        IDR residues / sequence length.
    """
    d = meta.predict_disorder_domains(seq, version=METAPREDICT_VERSION)
    return np.sum([len(x) for x in d.disordered_domains]) / len(seq)
