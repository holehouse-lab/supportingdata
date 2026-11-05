"""Example CAID3 proteins whose short IDRs metapredict V3 misses (Figure S10).

Each example is one full-length CAID3 Disorder-PDB protein drawn as three tracks
sharing a residue axis: the metapredict V3 score (with metapredict's own IDR calls
shaded), the AlphaFold2 pLDDT, and the Disorder-PDB reference labels with the short
IDRs that V3 misses outlined. The track drawing reuses ``ground_truth_panels`` (the
Figure 2D/E code) so the colours match Figure 2. Beside the tracks is the protein's
AlphaFold2 model, rendered by ``render_structures.py`` in the standard AlphaFold
pLDDT colours with the missed short IDRs drawn as thick, opaque tubes.

A short IDR is a maximal run of 10-49 residues labelled disordered in the
reference (CAID3 has none shorter than 10). Each is classed by what flanks it:

- ``loop``: ordered (PDB-resolved) residues on both sides;
- ``terminal tail``: it reaches the N- or C-terminus;
- ``next to unscored``: anything else, i.e. at least one neighbour is a residue
  CAID3 leaves unscored (neither DisProt-disordered nor resolved in the PDB).

V3 misses an IDR when fewer than half of its residues score at or above 0.5 *and*
fewer than half fall in a metapredict IDR call, so the trace and the shading agree.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import FigureBase
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Rectangle

from caid3_panels import DISORDER_THRESHOLD
from caid_eval import target_trace
from figtools import SHARED_DATA_DIR, shared_data_path
from ground_truth_panels import (
    C_DISORDER,
    C_IDR,
    C_NOINFO,
    C_ORDER,
    DIS_STATE_COLOR,
    METAPREDICT_COLOR,
    _draw_track,
    _idr_domains,
    _track_label,
    get_uniprot,
)

SHORT_IDR_MIN = 10  # residues
SHORT_IDR_MAX = 49  # residues
MISSED_FRACTION = 0.5  # an IDR is missed if V3 calls fewer than this fraction of it disordered
MAX_PROTEIN_LENGTH = 600  # residues; beyond this a 10-residue IDR is too thin to see

# Examples per context, in proportion to each context's share of the short-IDR
# residues V3 misses on CAID3 (loops 53%, next to unscored 26%, terminal tails 21%).
EXAMPLES_PER_CONTEXT = {"loop": 5, "next to unscored": 3, "terminal tail": 2}
PANEL_LETTERS = "ABCDEFGHIJ"

# CAID3 distributes AlphaFold-pLDDT as a disorder score of 1 - pLDDT/100
# (AlphaFold-disorder; Piovesan, Monzon & Tosatto, Protein Sci 31:e4466, 2022).
PLDDT_SCALE = 100.0
# pLDDT above 70 is AlphaFold's "confident" band (Jumper et al., Nature 596:583, 2021).
PLDDT_CONFIDENT = 70.0
# Standard AlphaFold DB pLDDT colour bands: (PyMOL colour name, lower bound, upper
# bound (exclusive; the top band's is above 100 so pLDDT = 100 is included), legend
# label, colour).
ALPHAFOLD_BANDS = [
    ("af_very_high", 90.0, PLDDT_SCALE + 1.0, "pLDDT > 90", "#0053D6"),
    ("af_confident", 70.0, 90.0, "pLDDT 70–90", "#65CBF3"),
    ("af_low", 50.0, 70.0, "pLDDT 50–70", "#FFDB13"),
    ("af_very_low", 0.0, 50.0, "pLDDT < 50", "#FF7D45"),
]
# Structure renders (and models.tsv) written by render_structures.py.
RENDER_DIR = SHARED_DATA_DIR / "alphafold_renders"
WHITE_LEVEL = 0.99  # render pixels at least this bright in every channel are background
CROP_MARGIN_PX = 10

DISPROT_CACHE = shared_data_path("pdb_cache", "disprot")
FIGSIZE = (6.5, 1.9)  # inches
TRACK_HEIGHT_RATIOS = [2.0, 1.6, 0.5]
TRACKS_TO_STRUCTURE_WIDTH = [3.2, 1.8]
TRACE_LW = 1.2
PLDDT_LW = 0.8
THRESHOLD_LW = 0.5
MISSED_OUTLINE_LW = 1.0
HEADER_FONTSIZE = 7.5
LABEL_FONTSIZE = 8
TICK_FONTSIZE = 7
LEGEND_FONTSIZE = 7
LEGEND_COLUMNS = 4  # legend above a single panel


@dataclass
class ShortIDR:
    """One short reference IDR: 0-based half-open span, flanking context, V3 coverage."""

    start: int
    end: int
    context: str
    fraction_scored: float  # fraction of residues with V3 score >= 0.5
    fraction_in_call: float  # fraction of residues inside a metapredict IDR call

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def missed(self) -> bool:
        return self.fraction_scored < MISSED_FRACTION and self.fraction_in_call < MISSED_FRACTION


@dataclass
class Example:
    """A selected protein and its short IDRs."""

    target: str
    context: str
    short_idrs: List[ShortIDR]

    @property
    def missed(self) -> List[ShortIDR]:
        return [idr for idr in self.short_idrs if idr.missed]


def idr_runs(labels: np.ndarray) -> List[Tuple[int, int]]:
    """Maximal runs of residues labelled disordered, as 0-based half-open spans.

    Parameters
    ----------
    labels : np.ndarray
        Reference labels (1 disordered, 0 ordered, -1 unscored).

    Returns
    -------
    list of (int, int)
        ``(start, end)`` of every IDR.
    """
    padded = np.r_[0, (labels == 1).astype(int), 0]
    edges = np.flatnonzero(np.diff(padded))
    return list(zip(edges[::2].tolist(), edges[1::2].tolist()))


def flank_context(labels: np.ndarray, start: int, end: int) -> str:
    """Class an IDR by its flanks: ``loop``, ``terminal tail`` or ``next to unscored``."""
    if start == 0 or end == labels.size:
        return "terminal tail"
    if labels[start - 1] == 0 and labels[end] == 0:
        return "loop"
    return "next to unscored"


def short_idrs(sequence: str, labels: np.ndarray, v3_scores: np.ndarray) -> List[ShortIDR]:
    """Every short IDR of one target, with how much of it metapredict V3 calls disordered.

    Parameters
    ----------
    sequence : str
        Protein sequence.
    labels : np.ndarray
        Reference labels (1 disordered, 0 ordered, -1 unscored).
    v3_scores : np.ndarray
        metapredict V3 score of every residue.

    Returns
    -------
    list of ShortIDR
        IDRs of ``SHORT_IDR_MIN``-``SHORT_IDR_MAX`` residues, N- to C-terminal.
    """
    assert v3_scores.size == labels.size == len(sequence), "scores, labels and sequence differ in length"
    assert np.isfinite(v3_scores).all(), "metapredict V3 must score every residue"
    in_call = np.zeros(labels.size, dtype=bool)
    for call_start, call_end in _idr_domains(sequence, v3_scores):
        in_call[call_start:call_end] = True
    found = []
    for start, end in idr_runs(labels):
        if not SHORT_IDR_MIN <= end - start <= SHORT_IDR_MAX:
            continue
        found.append(ShortIDR(
            start=start, end=end, context=flank_context(labels, start, end),
            fraction_scored=float(np.mean(v3_scores[start:end] >= DISORDER_THRESHOLD)),
            fraction_in_call=float(np.mean(in_call[start:end])),
        ))
    return found


def select_examples(reference: Dict[str, tuple], v3: Dict[str, Dict[int, float]]) -> List[Example]:
    """Pick the example proteins: the clearest cases of each context, in proportion.

    Each protein of at most ``MAX_PROTEIN_LENGTH`` residues with at least one missed
    short IDR is assigned the context holding most of its missed short-IDR residues.
    Within each context, proteins are ranked by missed short-IDR residues (then by
    number of missed short IDRs, then by ID) and the top ``EXAMPLES_PER_CONTEXT``
    are taken.

    Parameters
    ----------
    reference : dict
        Target -> ``(sequence, labels)``.
    v3 : dict
        metapredict V3 predictions, target -> ``{residue index: score}``.

    Returns
    -------
    list of Example
        Loops first, then next-to-unscored, then terminal tails.
    """
    candidates: Dict[str, List[Tuple[int, int, str, Example]]] = {c: [] for c in EXAMPLES_PER_CONTEXT}
    for target, (sequence, labels) in reference.items():
        if len(sequence) > MAX_PROTEIN_LENGTH:
            continue
        _positions, scores = target_trace(labels, v3[target])
        idrs = short_idrs(sequence, labels, scores)
        missed = [idr for idr in idrs if idr.missed]
        if not missed:
            continue
        residues_by_context = {c: sum(i.length for i in missed if i.context == c) for c in EXAMPLES_PER_CONTEXT}
        context = max(residues_by_context, key=lambda c: residues_by_context[c])
        n_missed_residues = sum(i.length for i in missed)
        candidates[context].append((n_missed_residues, len(missed), target, Example(target, context, idrs)))
    selected = []
    for context, n_examples in EXAMPLES_PER_CONTEXT.items():
        ranked = sorted(candidates[context], key=lambda row: (-row[0], -row[1], row[2]))
        assert len(ranked) >= n_examples, f"only {len(ranked)} candidates for context {context!r}"
        selected.extend(row[3] for row in ranked[:n_examples])
    return selected


def plddt_trace(labels: np.ndarray, alphafold_plddt: Dict[int, float]) -> np.ndarray:
    """AlphaFold2 pLDDT of every residue, recovered from the CAID3 AlphaFold-pLDDT score."""
    _positions, scores = target_trace(labels, alphafold_plddt)
    assert np.isfinite(scores).all(), "AlphaFold-pLDDT is missing residues for this target"
    return PLDDT_SCALE * (1.0 - scores)


def protein_name(target: str) -> str:
    """Protein name from the cached DisProt record."""
    with open(DISPROT_CACHE / f"{target}.json") as handle:
        return json.load(handle)["name"]


def describe_selection(examples: List[Example]) -> None:
    """Print each selected protein with its short IDRs and how many V3 misses."""
    for letter, example in zip(PANEL_LETTERS, examples):
        n_missed = len(example.missed)
        print(f"{letter}: {example.target} ({protein_name(example.target)}), {example.context}: "
              f"V3 misses {n_missed} of {len(example.short_idrs)} short IDRs")


def describe_example(example: Example, reference: Dict[str, tuple], alphafold_plddt: Dict[int, float]) -> None:
    """Print every short IDR of one example with V3 coverage and AlphaFold pLDDT.

    Parameters
    ----------
    example : Example
        The protein to describe.
    reference : dict
        Target -> ``(sequence, labels)``.
    alphafold_plddt : dict
        CAID3 AlphaFold-pLDDT scores of this target, ``{residue index: score}``.
    """
    sequence, labels = reference[example.target]
    plddt = plddt_trace(labels, alphafold_plddt)
    print(f"{example.target}  {protein_name(example.target)}  {len(sequence)} aa")
    print(f"{'residues':>10s} {'context':>17s} {'V3>=0.5':>8s} {'V3 IDR':>7s} {'pLDDT':>6s} {'missed':>6s}")
    for idr in example.short_idrs:
        span = f"{idr.start + 1}-{idr.end}"
        mean_plddt = float(np.mean(plddt[idr.start:idr.end]))
        print(f"{span:>10s} {idr.context:>17s} {idr.fraction_scored:8.0%} {idr.fraction_in_call:7.0%} "
              f"{mean_plddt:6.1f} {'yes' if idr.missed else 'no':>6s}")


def plddt_summary(reference: Dict[str, tuple], v3: Dict[str, Dict[int, float]],
                  alphafold_plddt: Dict[str, Dict[int, float]]) -> None:
    """Print how often the short IDRs V3 misses, and those it finds, are confidently folded by AlphaFold2.

    Covers every CAID3 target (no length limit). An IDR counts as confidently
    folded when its mean pLDDT is above ``PLDDT_CONFIDENT``.

    Parameters
    ----------
    reference : dict
        Target -> ``(sequence, labels)``.
    v3 : dict
        metapredict V3 predictions, target -> ``{residue index: score}``.
    alphafold_plddt : dict
        CAID3 AlphaFold-pLDDT predictions, target -> ``{residue index: score}``.
    """
    rows = []
    for target, (sequence, labels) in reference.items():
        _positions, scores = target_trace(labels, v3[target])
        plddt = plddt_trace(labels, alphafold_plddt[target])
        for idr in short_idrs(sequence, labels, scores):
            rows.append((idr.missed, idr.context, float(np.mean(plddt[idr.start:idr.end]))))
    print(f"{'':6s} {'context':17s} {'IDRs':>5s} {'median pLDDT':>13s} {'mean pLDDT > 70':>16s}")
    for missed in (True, False):
        for context in (None, *EXAMPLES_PER_CONTEXT):
            values = np.array([r[2] for r in rows if r[0] == missed and context in (None, r[1])])
            print(f"{'missed' if missed else 'found':6s} {context or 'all':17s} {values.size:5d} "
                  f"{np.median(values):13.1f} {np.mean(values > PLDDT_CONFIDENT):16.0%}")


def render_path(target: str) -> Path:
    """Path of the structure render of one target.

    Parameters
    ----------
    target : str
        DisProt ID of the CAID3 target, e.g. ``"DP04063"``.

    Returns
    -------
    pathlib.Path
        ``data/alphafold_renders/<target>.png``, written by ``render_structures.py``
        (the path is returned whether or not the render exists yet).
    """
    return RENDER_DIR / f"{target}.png"


def crop_to_content(image: np.ndarray) -> np.ndarray:
    """Trim the white background around a render, keeping a margin of ``CROP_MARGIN_PX``.

    Parameters
    ----------
    image : np.ndarray
        RGB(A) image with values in [0, 1], as read by ``matplotlib.image.imread``.

    Returns
    -------
    np.ndarray
        The cropped image.
    """
    content = np.any(image[..., :3] < WHITE_LEVEL, axis=2)
    rows = np.flatnonzero(content.any(axis=1))
    cols = np.flatnonzero(content.any(axis=0))
    assert rows.size and cols.size, "render is blank"
    top = max(int(rows[0]) - CROP_MARGIN_PX, 0)
    bottom = min(int(rows[-1]) + CROP_MARGIN_PX + 1, image.shape[0])
    left = max(int(cols[0]) - CROP_MARGIN_PX, 0)
    right = min(int(cols[-1]) + CROP_MARGIN_PX + 1, image.shape[1])
    return image[top:bottom, left:right]


def _legend_handles() -> list:
    """Legend entries for the three tracks and the structure render."""
    return [
        plt.Line2D([0], [0], color=METAPREDICT_COLOR, lw=TRACE_LW, label="metapredict V3"),
        Patch(facecolor=C_IDR, alpha=0.3, label="metapredict IDR call"),
        plt.Line2D([0], [0], color="black", lw=PLDDT_LW, label="AlphaFold2 pLDDT"),
        Patch(facecolor=C_DISORDER, label="disordered"),
        Patch(facecolor=C_ORDER, label="ordered"),
        Patch(facecolor=C_NOINFO, hatch="///", edgecolor="#b4b4b0", label="unscored"),
        Patch(facecolor="none", edgecolor="black", lw=MISSED_OUTLINE_LW,
              label="short IDR missed by V3\n(thick tube in the model)"),
        *[Patch(facecolor=colour, label=label) for _name, _lower, _upper, label, colour in ALPHAFOLD_BANDS],
    ]


def add_legend(fig: FigureBase, ncol: int = LEGEND_COLUMNS) -> None:
    """Draw the panel legend centred above ``fig``.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or matplotlib.figure.SubFigure
        The figure to put the legend above.
    ncol : int, optional
        Number of legend columns.
    """
    fig.legend(handles=_legend_handles(), loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=ncol,
               fontsize=LEGEND_FONTSIZE, frameon=False, columnspacing=1.2, handlelength=1.4)


def example_panel(example: Example, reference: Dict[str, tuple], v3: Dict[str, Dict[int, float]],
                  alphafold_plddt: Dict[str, Dict[int, float]], show_legend: bool = False,
                  fig: FigureBase | None = None) -> FigureBase:
    """Draw one example protein as three tracks beside its AlphaFold2 model and return the figure.

    Parameters
    ----------
    example : Example
        The protein to draw.
    reference : dict
        Target -> ``(sequence, labels)``.
    v3 : dict
        metapredict V3 predictions, target -> ``{residue index: score}``.
    alphafold_plddt : dict
        CAID3 AlphaFold-pLDDT predictions, target -> ``{residue index: score}``.
    show_legend : bool, optional
        Draw the legend above the panel.
    fig : matplotlib.figure.Figure or matplotlib.figure.SubFigure, optional
        Where to draw the panel, e.g. one subfigure of the combined Figure S10. A new
        figure of size ``FIGSIZE`` is made if not given.

    Returns
    -------
    matplotlib.figure.Figure or matplotlib.figure.SubFigure
        The panel (``fig`` if it was given). Its first axes is the top track.
    """
    sequence, labels = reference[example.target]
    length = len(sequence)
    positions, scores = target_trace(labels, v3[example.target])
    plddt = plddt_trace(labels, alphafold_plddt[example.target])

    render = render_path(example.target)
    assert render.exists(), f"{render} is missing: run `python render_structures.py` first"

    if fig is None:
        fig = plt.figure(figsize=FIGSIZE)
    grid = GridSpec(3, 2, figure=fig, height_ratios=TRACK_HEIGHT_RATIOS,
                    width_ratios=TRACKS_TO_STRUCTURE_WIDTH, hspace=0.35, wspace=0.05)
    ax_trace = fig.add_subplot(grid[0, 0])
    ax_plddt = fig.add_subplot(grid[1, 0], sharex=ax_trace)
    ax_labels = fig.add_subplot(grid[2, 0], sharex=ax_trace)
    ax_structure = fig.add_subplot(grid[:, 1])
    for ax in (ax_trace, ax_plddt, ax_labels):
        ax.tick_params(labelsize=TICK_FONTSIZE)

    # Track 1: metapredict V3 score with metapredict's IDR calls shaded (as in Figure 2D/E).
    for call_start, call_end in _idr_domains(sequence, scores):
        ax_trace.axvspan(call_start + 0.5, call_end + 0.5, color=C_IDR, alpha=0.16, lw=0, zorder=1)
    ax_trace.plot(positions, scores, color=METAPREDICT_COLOR, lw=TRACE_LW, zorder=4)
    ax_trace.axhline(DISORDER_THRESHOLD, ls="--", color="black", lw=THRESHOLD_LW, zorder=2)
    ax_trace.set_ylim(-0.02, 1.02)
    ax_trace.set_yticks([0, 0.5, 1])
    ax_trace.set_xlim(1, length)
    ax_trace.tick_params(labelbottom=False)
    _track_label(ax_trace, "metapredict\nV3 score")

    # Track 2: AlphaFold2 pLDDT, with the confident-band boundary.
    ax_plddt.plot(positions, plddt, color="black", lw=PLDDT_LW, zorder=4)
    ax_plddt.axhline(PLDDT_CONFIDENT, ls="--", color="black", lw=THRESHOLD_LW, zorder=2)
    ax_plddt.set_ylim(0, PLDDT_SCALE)
    ax_plddt.set_yticks([0, 50, 100])
    ax_plddt.tick_params(labelbottom=False)
    _track_label(ax_plddt, "AlphaFold2\npLDDT")

    # Track 3: Disorder-PDB labels, with the short IDRs V3 misses outlined.
    _draw_track(ax_labels, labels, DIS_STATE_COLOR)
    for idr in example.missed:
        ax_labels.add_patch(Rectangle((idr.start + 0.5, 0.0), idr.length, 1.0, fill=False,
                                      edgecolor="black", lw=MISSED_OUTLINE_LW, zorder=6, clip_on=False))
    _track_label(ax_labels, "Disorder-PDB")
    ax_labels.set_xlabel("Residue position", fontsize=LABEL_FONTSIZE)

    # AlphaFold2 model, pLDDT-coloured, missed short IDRs as thick tubes.
    # interpolation='none' embeds the render in the PDF at full resolution instead of
    # resampling it to the figure dpi.
    ax_structure.imshow(crop_to_content(mpimg.imread(render)), interpolation="none")
    ax_structure.axis("off")

    accession = get_uniprot(example.target)
    ident = f"{example.target} ({accession[0]})" if accession else example.target
    ax_trace.text(0.0, 1.04, f"{ident}  ·  {protein_name(example.target)}  ·  {length} aa",
                  transform=ax_trace.transAxes, fontsize=HEADER_FONTSIZE, ha="left", va="bottom")
    if show_legend:
        add_legend(fig)
    return fig
