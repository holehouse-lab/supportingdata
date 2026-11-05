"""Plotting helper for Figure S1: metapredict V1 disorder above AF2 pLDDT."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from figtools import V1_COLOR

# metapredict V1 disorder threshold (Emenecker et al. 2021, Biophys. J.)
V1_DISORDER_THRESHOLD = 0.42
# AF2 pLDDT below 70 is "low" or "very low" confidence (Jumper et al. 2021, Nature)
PLDDT_LOW_CONFIDENCE = 70
PLDDT_COLOR = "black"


def shade_runs(ax: plt.Axes, mask: np.ndarray, color: str, alpha: float) -> None:
    """Shade each contiguous run of residues where ``mask`` is True.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to shade.
    mask : numpy.ndarray of bool
        One entry per residue (residue i is at x = i, 1-indexed).
    color : str
        Shading colour.
    alpha : float
        Shading opacity.
    """
    padded = np.concatenate([[False], mask, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(int)))
    for start, end in zip(edges[::2], edges[1::2]):
        # run covers residues start+1 .. end (1-indexed)
        ax.axvspan(start + 0.5, end + 0.5, color=color, alpha=alpha, linewidth=0, zorder=0)


def plot_disorder_vs_plddt(disorder: np.ndarray, plddt: np.ndarray, label: str) -> plt.Figure:
    """Two-row panel: V1 disorder (top) and AF2 pLDDT (bottom) for one protein.

    Residues above the V1 disorder threshold (top) and below pLDDT 70 (bottom)
    are shaded, so the two methods' calls can be compared directly.

    Parameters
    ----------
    disorder : numpy.ndarray
        metapredict V1 per-residue disorder scores (0-1).
    plddt : numpy.ndarray
        AF2 per-residue pLDDT (0-100).
    label : str
        Protein identifier written above the top axes, right-aligned.

    Returns
    -------
    matplotlib.figure.Figure
        The panel.
    """
    assert len(disorder) == len(plddt), "disorder and pLDDT must cover the same residues"
    residue = np.arange(1, len(disorder) + 1)

    fig, (ax_d, ax_p) = plt.subplots(nrows=2, figsize=(3.25, 3.5), dpi=150, sharex=True)

    shade_runs(ax_d, disorder >= V1_DISORDER_THRESHOLD, V1_COLOR, alpha=0.2)
    ax_d.plot(residue, disorder, linewidth=2, color=V1_COLOR)
    ax_d.axhline(V1_DISORDER_THRESHOLD, linestyle="dashed", color="black")
    ax_d.set_ylim([0, 1])
    ax_d.set_yticks([0, 0.5, 1], [0, 0.5, 1])
    ax_d.set_ylabel("Disorder score")
    ax_d.text(1.0, 1.02, label, transform=ax_d.transAxes, ha="right", va="bottom", fontsize=9)

    shade_runs(ax_p, plddt <= PLDDT_LOW_CONFIDENCE, PLDDT_COLOR, alpha=0.12)
    ax_p.plot(residue, plddt, linewidth=2, color=PLDDT_COLOR)
    ax_p.axhline(PLDDT_LOW_CONFIDENCE, linestyle="dashed", color="black")
    ax_p.set_ylim([0, 100])
    ax_p.set_ylabel("AF2 pLDDT")
    ax_p.set_xlabel("Residue")
    ax_p.set_xlim([1, len(disorder)])

    plt.tight_layout()
    return fig
