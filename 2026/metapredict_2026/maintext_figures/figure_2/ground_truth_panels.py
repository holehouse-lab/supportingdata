"""Best/worst-predicted protein panels for Figure 2D/E (four tracks per protein).

Adapted from ``CAID3_comparison/figures/ground_truth_traces.py``. Each selected
protein is drawn as a block of four tracks sharing a residue axis:

1. **metapredict V3** score, with metapredict's own IDR calls (domain
   decomposition applied to these scores) shaded;
2. **other CAID3 predictors** (PUNCH2, SPOT-Disorder2, AlphaFold-rsa);
3. **PDB** - the strongest structural evidence for each residue: structured as a
   monomer > as a homo-oligomer > only in a hetero-complex, or no electron density
   / not in any structure;
4. **Disorder-PDB** - the CAID3 reference labels.

Changes from the original: metapredict V3 comes from the 3.1.0 predictions; the
comparison predictors use the fixed manuscript colours; the PDB-coverage cache is
read from the shared data directory (built by ``data/pdb_cache_pipeline``), with
the structural-state codes and UniProt lookup inlined instead of importing the
network-fetching pipeline; the figure suptitle is dropped and the composite is
returned rather than saved; IDR-domain failures raise instead of being skipped.

Restyled in November 2026 to match the other manuscript figures: the composite is
drawn at half page width with the text sizes of Figures 1, S10 and S11; the PDB track
is drawn in greys (so it shares no colour with a predictor line); the PDB codes are
printed by :func:`print_pdb_codes` instead of being written under each block; track
labels are drawn on the first block only; and the NMR Cα RMSF overlay is optional
(off in Figure 2).
"""

from __future__ import annotations

import json
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from metapredict.backend.domain_definition import get_domains
from sklearn.metrics import roc_auc_score

from caid_eval import target_trace
from figtools import (
    CAID3_PREDICTOR_COLORS,
    DISORDER_PDB_COLORS,
    HATCH_EDGE_COLOR,
    V3_COLOR,
    shared_data_path,
)

METAPREDICT_COLOR = V3_COLOR

# Per-residue structural states written by data/pdb_cache_pipeline/structure/pdb_coverage.py
ST_UNKNOWN = 0          # not covered by any PDB construct
ST_NODENSITY = 1        # in a construct but unresolved anywhere (missing density)
ST_STRUCTURED_MONO = 2  # resolved as a monomer (lone chain) in at least one structure
ST_COMPLEX = 3          # resolved only in a hetero-complex (never as monomer or homo)
ST_STRUCTURED_HOMO = 4  # resolved as a homo-oligomer, but never as a monomer
ST_NODATA = -1          # no PDB data retrievable

PDB_CACHE = shared_data_path("pdb_cache")
COVERAGE_JSON = PDB_CACHE / "coverage_states.json"
COVERAGE_PDBS_JSON = PDB_CACHE / "coverage_pdbs.json"
RMSF_JSON = PDB_CACHE / "rmsf.json"
DISPROT_CACHE = PDB_CACHE / "disprot"

# NMR RMSF gradient: resolved residues from an NMR ensemble are coloured by their
# per-residue Cα RMSF on a fixed 0–8 Å blue→red scale (values >8 saturate red), so
# the colour is comparable across every protein.
RMSF_CAP = 8.0
RMSF_CMAP = LinearSegmentedColormap.from_list(
    "rmsf", ["#2f62ab", "#7ba7d6", "#f0efe6", "#e08a5a", "#c0392b"])
RMSF_CMAP.set_bad(alpha=0.0)
RMSF_NORM = Normalize(vmin=0.0, vmax=RMSF_CAP)

# Disorder-label track colours (Disorder-PDB), shared with Figures 2C, S4 and S10.
C_DISORDER = DISORDER_PDB_COLORS[1]   # orange
C_ORDER = DISORDER_PDB_COLORS[0]      # green
C_NOINFO = DISORDER_PDB_COLORS[-1]    # light grey (hatched) — unknown / excluded / not observed

# PDB-coverage track colours. Resolved residues are greys, darkest for the strongest
# evidence of folding on its own (a lone monomer), so the track shares no hue with
# the predictor lines or the Disorder-PDB bar; missing density, a disorder cue, is a
# hatched pale version of the Disorder-PDB "disordered" orange.
C_STRUCTURED_MONO = "#3c3c3c"  # dark grey  — resolved as a monomer in isolation
C_STRUCTURED_HOMO = "#7d7d7d"  # mid grey   — resolved in isolation as a homo-oligomer
C_COMPLEX = "#b4b4b4"          # light grey — resolved only in a hetero-complex
C_NODENSITY = "#f4c4a0"        # pale orange (hatched) — in a construct, missing density

# metapredict IDR-domain shading on the trace (a pale wash of the metapredict V3
# colour, since it marks metapredict's own predicted IDR domains).
C_IDR = METAPREDICT_COLOR
IDR_ALPHA = 0.16

# Additional predictors shown as a comparison trace under metapredict (top CAID3
# Disorder-PDB methods). Name (== .caid filename) -> line colour.
OTHER_PREDICTORS = [(name, CAID3_PREDICTOR_COLORS[name])
                    for name in ("PUNCH2", "SPOT-Disorder2", "AlphaFold-rsa")]

# Composite layout (inches). Each protein is a block of four tracks sharing a residue
# axis. The track width is chosen so D (with track labels) and E (without) fit side by
# side across a 6.4 in page.
AXES_WIDTH_IN = 2.5
TRACK_HEIGHTS_IN = [0.36, 0.26, 0.07, 0.07]  # V3 trace, other predictors, PDB, Disorder-PDB
TRACK_GAPS_IN = [0.1, 0.07, 0.05]            # between consecutive tracks (room for the 0/1 labels)
BLOCK_GAP_IN = 0.34                          # x tick labels + the next protein's header
ROWS_PER_BLOCK = 2 * len(TRACK_HEIGHTS_IN)   # tracks, the gaps between them, the block gap
TRACK_LABELS = ["metapredict\nV3", "other\nCAID3", "PDB", "Disorder-PDB"]
V3_LW = 1.2
OTHER_LW = 0.8
THRESHOLD_LW = 0.5

# Text sizes (points), matching Figures 1, S10 and S11.
TICK_FONTSIZE = 7
LABEL_FONTSIZE = 8
HEADER_FONTSIZE = 7
LEGEND_FONTSIZE = 6.5
LEGEND_COLUMNS = 3

DIS_STATE_COLOR = {1: C_DISORDER, 0: C_ORDER, -1: C_NOINFO}
COV_STATE_COLOR = {
    ST_STRUCTURED_MONO: C_STRUCTURED_MONO,
    ST_STRUCTURED_HOMO: C_STRUCTURED_HOMO,
    ST_COMPLEX: C_COMPLEX,
    ST_NODENSITY: C_NODENSITY,
    ST_UNKNOWN: C_NOINFO,
}

# Resolved structural states (any electron density), used for the RMSF overlay.
_RESOLVED_STATES = (ST_STRUCTURED_MONO, ST_STRUCTURED_HOMO, ST_COMPLEX)

# Legend / printout names of the resolved states, keyed like coverage_pdbs.json.
REGION_KINDS = [(ST_STRUCTURED_MONO, "mono", "monomer"),
                (ST_STRUCTURED_HOMO, "homo", "homo-oligomer"),
                (ST_COMPLEX, "cplx", "hetero-complex only")]


def get_uniprot(disprot_id: str) -> Optional[Tuple[str, str]]:
    """Return ``(uniprot_acc, sequence)`` for a DisProt target from the cached API response."""
    with open(DISPROT_CACHE / f"{disprot_id}.json") as handle:
        doc = json.load(handle)
    if not doc or "acc" not in doc:
        return None
    return doc["acc"], doc.get("sequence", "")


def _runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return (start, end) inclusive index runs where ``mask`` is True."""
    if not mask.any():
        return []
    idx = np.flatnonzero(mask)
    runs = []
    start = idx[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if b != a + 1:
            runs.append((start, a))
            start = b
    runs.append((start, idx[-1]))
    return runs


def _draw_track(
    ax, labels: np.ndarray, color_map: Dict[int, str], hatch_states=(-1,)
) -> None:
    """Fill a thin annotation track with coloured blocks per state value."""
    for state, color in color_map.items():
        hatch = "///" if state in hatch_states else None
        for s, e in _runs(labels == state):
            ax.axvspan(
                s + 0.5, e + 1.5, facecolor=color, lw=0,
                hatch=hatch, edgecolor=HATCH_EDGE_COLOR if hatch else "none",
            )
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _region_codes(state: np.ndarray, state_value: int, observing_runs: List) -> List[dict]:
    """Return the resolved runs of ``state_value`` with their PDB code(s).

    Each element is ``{"start", "end", "codes"}`` (0-based inclusive indices, PDB
    codes upper-cased and sorted). Runs with no associated PDB code are dropped.
    """
    if not observing_runs:
        return []
    L = state.size
    per_res: List[set] = [set() for _ in range(L)]
    for start, end, ids in observing_runs:  # 1-based inclusive
        for pos in range(start - 1, min(end, L)):
            per_res[pos].update(ids)
    regions = []
    for s, e in _runs(state == state_value):
        ids: set = set()
        for pos in range(s, e + 1):
            ids |= per_res[pos]
        if not ids:
            continue
        regions.append({"start": s, "end": e, "codes": sorted(i.upper() for i in ids)})
    return regions


def pdb_code_regions(cov_states: np.ndarray, cov_pdbs: Optional[dict]) -> List[dict]:
    """List each resolved region of one protein's PDB track with the PDB codes behind it.

    Parameters
    ----------
    cov_states : numpy.ndarray
        Per-residue structural states (``ST_*`` codes) from ``coverage_states.json``.
    cov_pdbs : dict or None
        The protein's entry in ``coverage_pdbs.json``: ``{"mono" | "homo" | "cplx":
        [[first, last, [codes]], ...]}`` with 1-based inclusive residue ranges.

    Returns
    -------
    list of dict
        One ``{"start", "end", "kind", "codes"}`` per resolved run, left to right;
        ``start``/``end`` are 1-based inclusive residues and ``kind`` is the legend
        name of the structural state.
    """
    cov_pdbs = cov_pdbs or {}
    regions = []
    for state_value, key, kind in REGION_KINDS:
        for r in _region_codes(cov_states, state_value, cov_pdbs.get(key, [])):
            regions.append({"start": r["start"] + 1, "end": r["end"] + 1, "kind": kind,
                            "codes": r["codes"]})
    return sorted(regions, key=lambda r: r["start"])


def print_pdb_codes(selection, coverage: Dict[str, np.ndarray],
                    coverage_pdbs: Dict[str, dict]) -> None:
    """Print the PDB codes behind every resolved block drawn on the PDB track.

    The codes used to be written under each protein in the figure; they are printed
    here instead (for the caption or Methods) to keep the panel free of text.

    Parameters
    ----------
    selection : list of tuple
        ``(target_id, auc)`` pairs, as returned by :func:`select_best_worst`.
    coverage : dict
        Target -> per-residue structural states, from :func:`load_coverage`.
    coverage_pdbs : dict
        Target -> PDB codes per state, from :func:`load_coverage`.
    """
    for tid, _auc in selection:
        states = coverage.get(tid)
        if states is None:
            print(f"{tid}: no PDB data retrieved")
            continue
        regions = pdb_code_regions(states, coverage_pdbs.get(tid))
        listing = "; ".join(f"{r['start']}-{r['end']} {r['kind']}: {', '.join(r['codes'])}"
                            for r in regions)
        print(f"{tid}: {listing or 'no resolved residues'}")


def _gradient_block(ax, rmsf: Optional[np.ndarray], resolved: np.ndarray,
                    L: int) -> bool:
    """Recolour resolved NMR residues by per-residue Cα RMSF (blue→red, 0–8 Å).

    Overpaints the categorical block colour on residues that are both resolved
    (structured/complex) and covered by an NMR ensemble; other residues are left
    transparent so their categorical colour shows through. Returns True if any
    residue was recoloured.
    """
    if rmsf is None:
        return False
    show = np.isfinite(rmsf) & resolved
    if not show.any():
        return False
    grid = np.full(L, np.nan)
    grid[show] = rmsf[show]
    ax.imshow(
        np.ma.masked_invalid(grid.reshape(1, L)), aspect="auto",
        extent=(0.5, L + 0.5, 0, 1), cmap=RMSF_CMAP, norm=RMSF_NORM,
        interpolation="nearest", zorder=5,
    )
    return True


def _idr_domains(seq: str, scores: np.ndarray) -> List[Tuple[int, int]]:
    """metapredict IDR domain boundaries from these scores.

    Applies metapredict's V3 domain-decomposition (savgol smoothing, 0.5
    threshold, 12-residue minimum IDR, 50-residue minimum folded domain,
    10-residue gap closure) directly to the supplied per-residue scores, so the
    shaded IDRs are consistent with the plotted trace. ``use_python=True`` avoids
    the optional Cython path. Returns ``[start, end)`` spans (0-based).
    """
    vals = np.where(np.isnan(scores), 0.0, scores).astype(np.float64)
    _smoothed, idrs, _fds = get_domains(
        seq, vals, disorder_threshold=0.5, minimum_IDR_size=12,
        minimum_folded_domain=50, gap_closure=10, use_python=True,
    )
    return [(idr[0], idr[1]) for idr in idrs]


def _track_label(ax, text: str, fontsize: float = 8) -> None:
    """Write a track's name horizontally to the left of it (as in Figure S10).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The track.
    text : str
        Track name.
    fontsize : float, optional
        Text size in points (Figure S10 uses the default of 8).
    """
    ax.set_ylabel(text, rotation=0, ha="right", va="center", fontsize=fontsize)


def _score_track(ax) -> None:
    """Shared look of a 0-1 score track: dashed 0.5 threshold, ticks at 0 and 1."""
    ax.axhline(0.5, ls="--", color="black", lw=THRESHOLD_LW, zorder=2)
    ax.set_ylim(-0.02, 1.02)
    ax.set_yticks([0, 1])
    ax.tick_params(labelsize=TICK_FONTSIZE, labelbottom=False)


def _render_tracks(
    axes,
    seq: str,
    pdb_labels: np.ndarray,
    cov_states: Optional[np.ndarray],
    rmsf: Optional[np.ndarray],
    mp_pos: np.ndarray,
    mp_score: np.ndarray,
    others: Dict[str, np.ndarray],
    label_tracks: bool,
) -> bool:
    """Draw one protein's four tracks onto ``axes``.

    Parameters
    ----------
    axes : tuple of matplotlib.axes.Axes
        ``(ax_trace, ax_other, ax_cov, ax_pdb)``, sharing the residue axis.
    seq : str
        Protein sequence.
    pdb_labels : numpy.ndarray
        Disorder-PDB labels (1 / 0 / -1) per residue.
    cov_states : numpy.ndarray or None
        Per-residue structural states (``ST_*``); ``None`` if no PDB data.
    rmsf : numpy.ndarray or None
        Per-residue NMR Cα RMSF (Å); ``None`` draws no RMSF overlay.
    mp_pos, mp_score : numpy.ndarray
        Residue positions and metapredict V3 scores.
    others : dict
        Predictor name -> per-residue scores for the comparison track.
    label_tracks : bool
        Write the track names to the left (done for the first protein only).

    Returns
    -------
    bool
        True if any residue was recoloured by the NMR RMSF gradient.
    """
    ax_trace, ax_other, ax_cov, ax_pdb = axes
    L = len(seq)

    # Track 1: metapredict trace + red IDR domain shading. IDR spans are 0-based
    # half-open [start, end); residue (start+1..end) maps to x-span (start+0.5, end+0.5).
    for s, e in _idr_domains(seq, mp_score):
        ax_trace.axvspan(s + 0.5, e + 0.5, color=C_IDR, alpha=IDR_ALPHA, lw=0, zorder=1)
    ax_trace.plot(mp_pos, mp_score, color=METAPREDICT_COLOR, lw=V3_LW, zorder=4)
    _score_track(ax_trace)
    ax_trace.set_xlim(1, L)

    # Track 2: comparison predictors (top CAID3 Disorder-PDB methods).
    for name, color in OTHER_PREDICTORS:
        s = others.get(name)
        if s is None or not np.isfinite(s).any():
            continue
        ax_other.plot(np.arange(1, L + 1), s, color=color, lw=OTHER_LW, zorder=4)
    _score_track(ax_other)

    # Track 3: PDB structural coverage.
    if cov_states is not None and not (cov_states == ST_NODATA).all():
        _draw_track(ax_cov, cov_states, COV_STATE_COLOR,
                    hatch_states=(ST_UNKNOWN, ST_NODENSITY))
        resolved = np.isin(cov_states, _RESOLVED_STATES)
        has_rmsf = _gradient_block(ax_cov, rmsf, resolved, L)
    else:
        has_rmsf = False
        ax_cov.text(0.5, 0.5, "no PDB data", transform=ax_cov.transAxes,
                    ha="center", va="center", fontsize=TICK_FONTSIZE)
        ax_cov.set_yticks([])
        for spine in ax_cov.spines.values():
            spine.set_visible(False)
    ax_cov.tick_params(bottom=False, labelbottom=False)

    # Track 4: Disorder-PDB reference.
    _draw_track(ax_pdb, pdb_labels, DIS_STATE_COLOR)
    ax_pdb.tick_params(labelsize=TICK_FONTSIZE)

    if label_tracks:
        for ax, text in zip(axes, TRACK_LABELS):
            _track_label(ax, text, fontsize=TICK_FONTSIZE)
    return has_rmsf


def _legend_handles() -> list:
    """Legend entries for what Figure 2A-C does not already key: IDR calls and the PDB track.

    The predictor colours are keyed in Figure 2B and the Disorder-PDB colours in
    Figure 2C.
    """
    return [
        Patch(facecolor=C_IDR, alpha=IDR_ALPHA * 2, label="V3 IDR call"),
        Patch(facecolor=C_STRUCTURED_MONO, label="monomer"),
        Patch(facecolor=C_STRUCTURED_HOMO, label="homo-oligomer"),
        Patch(facecolor=C_COMPLEX, label="hetero-complex"),
        Patch(facecolor=C_NODENSITY, hatch="///", edgecolor=HATCH_EDGE_COLOR,
              label="missing density"),
        Patch(facecolor=C_NOINFO, hatch="///", edgecolor=HATCH_EDGE_COLOR,
              label="no structure"),
    ]


def _per_protein_auc(pdb_labels: np.ndarray, pred_scores: Dict[int, float]):
    """Per-protein ROC-AUC of a predictor over the scored residues of one target.

    Returns ``(auc, n_disordered, n_ordered)``; ``auc`` is ``None`` when the target
    has only one class among its scored, covered residues (AUC undefined).
    """
    y_true, y_score = [], []
    for pos, label in enumerate(pdb_labels, start=1):
        if label < 0:
            continue
        score = pred_scores.get(pos)
        if score is None:
            continue
        y_true.append(int(label))
        y_score.append(score)
    y_true = np.asarray(y_true)
    n_dis = int((y_true == 1).sum())
    n_ord = int((y_true == 0).sum())
    if n_dis == 0 or n_ord == 0:
        return None, n_dis, n_ord
    return float(roc_auc_score(y_true, np.asarray(y_score))), n_dis, n_ord


def select_best_worst(
    reference: Dict[str, tuple], mp: Dict[str, Dict[int, float]], n: int = 5,
    min_per_class: int = 15, min_len: int = 40, max_len: int = 500,
):
    """Rank targets by metapredict V3 per-protein AUC and return ``(best, worst)``.

    Only targets with at least ``min_per_class`` disordered *and* ordered scored
    residues, within ``[min_len, max_len]`` residues, are considered, so each panel
    shows a substantive mix of order and disorder. Each returned element is
    ``(target_id, auc)``. Best is highest AUC first; worst is lowest AUC first.
    Ties break toward more disordered residues (the more illustrative case).
    """
    rows = []
    for tid, (seq, lab) in reference.items():
        if not (min_len <= len(seq) <= max_len):
            continue
        auc, n_dis, n_ord = _per_protein_auc(lab, mp.get(tid, {}))
        if auc is None or n_dis < min_per_class or n_ord < min_per_class:
            continue
        rows.append((tid, auc, n_dis))
    worst = sorted(rows, key=lambda r: (r[1], -r[2]))[:n]
    best = sorted(rows, key=lambda r: (-r[1], -r[2]))[:n]
    return [(t, a) for t, a, _ in best], [(t, a) for t, a, _ in worst]


def composite_height(n_blocks: int) -> float:
    """Height in inches of the track area of a :func:`make_composite` figure.

    Parameters
    ----------
    n_blocks : int
        Number of proteins.

    Returns
    -------
    float
        Height in inches; the legend, headers and tick labels sit outside it and
        are kept by the tight bounding box used when saving.
    """
    block = sum(TRACK_HEIGHTS_IN) + sum(TRACK_GAPS_IN)
    return n_blocks * block + (n_blocks - 1) * BLOCK_GAP_IN


def _block_height_ratios(n_blocks: int) -> List[float]:
    """GridSpec height ratios (inches): per block, each track then the gap after it."""
    ratios: List[float] = []
    for i in range(n_blocks):
        for height, gap in zip(TRACK_HEIGHTS_IN, TRACK_GAPS_IN + [BLOCK_GAP_IN]):
            ratios.append(height)
            ratios.append(gap)
    return ratios[:-1]  # no gap after the last block


def make_composite(
    selection,
    reference: Dict[str, tuple], coverage: Dict[str, np.ndarray],
    mp: Dict[str, Dict[int, float]], other_preds: Dict[str, Dict[str, Dict[int, float]]],
    rmsf: Optional[Dict[str, np.ndarray]] = None, legend: bool = True,
    track_labels: bool = True,
) -> plt.Figure:
    """Draw the selected proteins as stacked four-track blocks and return the figure.

    Each block is headed by the protein's DisProt ID and metapredict V3 per-protein
    AUC; the track names are written beside the first block only, since every block
    has the same four tracks.

    Parameters
    ----------
    selection : list of tuple
        ``(target_id, auc)`` pairs, top to bottom (from :func:`select_best_worst`).
    reference : dict
        Target -> ``(sequence, labels)``.
    coverage : dict
        Target -> per-residue structural states, from :func:`load_coverage`.
    mp : dict
        metapredict V3 predictions, target -> ``{residue index: score}``.
    other_preds : dict
        Comparison predictor name -> its predictions, as ``mp``.
    rmsf : dict, optional
        Target -> per-residue NMR Cα RMSF (Å), from :func:`load_rmsf`. When given,
        resolved NMR residues on the PDB track are recoloured on a 0-8 Å blue-red
        scale and a colourbar is added. Figure 2 leaves it out.
    legend : bool, optional
        Draw the legend (IDR calls and the PDB-track states) above the panel.
        Figure 2 draws it on D only, since E has the same tracks.
    track_labels : bool, optional
        Name the four tracks beside the first block. Figure 2 does this on D only,
        since E sits to its right with the same tracks.

    Returns
    -------
    matplotlib.figure.Figure
        The composite; its tracks are ``AXES_WIDTH_IN`` wide.
    """
    n = len(selection)
    height_in = composite_height(n)
    fig = plt.figure(figsize=(AXES_WIDTH_IN, height_in))
    grid = GridSpec(ROWS_PER_BLOCK * n - 1, 1, figure=fig,
                    height_ratios=_block_height_ratios(n), hspace=0,
                    left=0.0, right=1.0, top=1.0, bottom=0.0)

    any_rmsf = False
    for i, (tid, auc) in enumerate(selection):
        seq, pdb_labels = reference[tid]
        first_row = ROWS_PER_BLOCK * i
        ax_trace = fig.add_subplot(grid[first_row, 0])
        ax_other = fig.add_subplot(grid[first_row + 2, 0], sharex=ax_trace)
        ax_cov = fig.add_subplot(grid[first_row + 4, 0], sharex=ax_trace)
        ax_pdb = fig.add_subplot(grid[first_row + 6, 0], sharex=ax_trace)

        pos, score = target_trace(pdb_labels, mp.get(tid, {}))
        others = {name: target_trace(pdb_labels, preds.get(tid, {}))[1]
                  for name, preds in other_preds.items()}
        has = _render_tracks(
            (ax_trace, ax_other, ax_cov, ax_pdb), seq, pdb_labels, coverage.get(tid),
            rmsf.get(tid) if rmsf else None, pos, score, others,
            label_tracks=track_labels and i == 0,
        )
        any_rmsf = any_rmsf or has
        ax_trace.text(0.0, 1.06, f"{tid}  ·  AUC = {auc:.2f}", transform=ax_trace.transAxes,
                      fontsize=HEADER_FONTSIZE, ha="left", va="bottom")
    ax_pdb.set_xlabel("Residue", fontsize=LABEL_FONTSIZE)

    if legend:
        # above the first protein's header line, right-aligned with the tracks so any
        # extra width runs over the (empty) track-label column
        header_top_in = 0.06 * TRACK_HEIGHTS_IN[0] + 0.14
        fig.legend(handles=_legend_handles(), loc="lower right", ncol=LEGEND_COLUMNS,
                   bbox_to_anchor=(1.0, 1.0 + header_top_in / height_in), frameon=False,
                   fontsize=LEGEND_FONTSIZE, columnspacing=1.0, handlelength=1.4,
                   borderaxespad=0, borderpad=0)

    if any_rmsf:
        cax = fig.add_axes((1.03, 0.35, 0.02, 0.30))
        sm = ScalarMappable(norm=RMSF_NORM, cmap=RMSF_CMAP)
        sm.set_array([])
        cb = fig.colorbar(cax=cax, mappable=sm)
        cb.set_ticks([0, 4, 8])
        cb.set_ticklabels(["0", "4", "≥8"])
        cax.tick_params(labelsize=TICK_FONTSIZE, length=2, pad=1)
        cb.set_label("NMR Cα RMSF (Å)", fontsize=LEGEND_FONTSIZE, labelpad=3)

    return fig


def load_coverage() -> Tuple[Dict[str, np.ndarray], Dict[str, List]]:
    """Load the per-residue PDB structural states and the PDB codes behind them.

    Returns
    -------
    tuple of dict
        ``(states, pdbs)``: target -> ``ST_*`` state per residue, and target ->
        ``{"mono" | "homo" | "cplx": [[first, last, [codes]], ...]}``.
    """
    with open(COVERAGE_JSON) as handle:
        states = {k: np.asarray(v, dtype=np.int8) for k, v in json.load(handle).items()}
    pdbs = {}
    with open(COVERAGE_PDBS_JSON) as handle:
        pdbs = json.load(handle)
    return states, pdbs


def load_rmsf() -> Dict[str, np.ndarray]:
    with open(RMSF_JSON) as handle:
        raw = json.load(handle)
    return {k: np.array([np.nan if v is None else v for v in vals], dtype=np.float64)
            for k, vals in raw.items()}
