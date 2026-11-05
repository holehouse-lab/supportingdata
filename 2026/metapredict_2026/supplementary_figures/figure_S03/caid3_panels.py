"""CAID3 Disorder-PDB loading and panel-drawing helpers.

Adapted from ``CAID3_comparison/figures/ranking_curves.py`` (the original CAID3
re-analysis). The scoring is unchanged; the differences are:

* metapredict V1/V2/V3 are scored from metapredict 3.1.0 predictions
  (``data/metapredict_caid`` in the shared data directory) instead of the
  metapredict files distributed by CAID3, which were made with earlier versions;
* every function draws onto an axes it is given, so each manuscript panel is its
  own figure, and titles/suptitles are dropped (identifiers go in-axis or on the
  y-axis);
* the top CAID3 methods keep one fixed colour in every panel
  (``figtools.CAID3_PREDICTOR_COLORS``) and metapredict uses the manuscript colours.

An identical copy of this file sits in every figure directory that uses it.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from caid_eval import (
    PredictorResult,
    collect_pooled,
    evaluate_predictor,
    list_predictor_files,
    parse_caid,
    parse_reference,
    pr_curve_xy,
    roc_curve_xy,
    target_trace,
)
from figtools import (
    CAID3_PREDICTOR_COLORS,
    DISORDER_PDB_COLORS,
    FIELD_COLOR,
    HATCH_EDGE_COLOR,
    VERSION_COLORS,
    metapredict_caid_path,
    shared_data_path,
)
from predictor_classes import is_disorder_predictor

DATASET_NAME = "Disorder-PDB"
# CAID3 Disorder-PDB reference set, as path components below the shared data directory.
REFERENCE_FASTA_PARTS = ("caid_references", "caid3_disorder_pdb.fasta")

# CAID3 file names of the metapredict versions -> manuscript labels.
METAPREDICT_NAMES = {
    "Metapredict-v1": "metapredict V1",
    "Metapredict-v2": "metapredict V2",
    "Metapredict-v3": "metapredict V3",
}
METAPREDICT_COLORS = {name: VERSION_COLORS[label.split()[-1]]
                      for name, label in METAPREDICT_NAMES.items()}

N_TOP = 5  # top-N CAID3 predictors (by AUC) drawn alongside metapredict

# Colours for any non-metapredict predictor outside CAID3_PREDICTOR_COLORS.
EXTRA_COLORS = ["#8a8a2a", "#d1189a", "#5b5b5b"]  # olive, magenta, dark grey

# Light grey for the unlabelled field of other predictors behind the highlighted
# curves (deliberately de-emphasised).
FIELD_CURVE_COLOR = "#d2d2d0"
FIELD_LABEL = "other CAID3 predictors"  # legend label for the grey field, as in Figure S11

DISORDER_THRESHOLD = 0.5

# Text sizes (points) for the curve and trace panels, matching Figures 1, S10 and S11.
TICK_FONTSIZE = 7
LABEL_FONTSIZE = 8
LEGEND_FONTSIZE = 6.5

# Legend labels for the Disorder-PDB bar, keyed like DISORDER_PDB_COLORS.
DISORDER_PDB_LABELS = {1: "disordered", 0: "ordered", -1: "unscored"}

# Per-residue trace layout (inches): each protein is a trace with its Disorder-PDB
# bar directly below it.
TRACE_HEIGHT_IN = 0.68
LABEL_BAR_HEIGHT_IN = 0.07
TRACE_TO_BAR_GAP_IN = 0.03
BLOCK_GAP_IN = 0.36  # room for the x tick labels and the next protein's ID
HEADER_HEIGHT_IN = 0.14  # height of the protein ID line above a trace
YLABEL_OFFSET_IN = 0.3  # shared y label sits this far left of the axes, clear of tick labels
LEGEND_COLUMNS = 6  # predictor + Disorder-PDB legend above a trace panel


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def display_name(name: str) -> str:
    """Manuscript label for a predictor (metapredict versions are relabelled)."""
    return METAPREDICT_NAMES.get(name, name)


def load_reference() -> Dict[str, Tuple[str, np.ndarray]]:
    """Load the CAID3 Disorder-PDB reference set (319 targets).

    The file is ``data/caid_references/caid3_disorder_pdb.fasta``; the path is
    resolved when this is called, so figures that never load the reference (S9)
    do not need it.

    Returns
    -------
    dict
        DisProt ID -> ``(sequence, labels)``, where ``labels`` is an integer array
        with 1 (disordered) or 0 (ordered) for scored residues and -1 for excluded
        residues (see ``caid_eval.parse_reference``).
    """
    return parse_reference(str(shared_data_path(*REFERENCE_FASTA_PARTS)))


def load_predictions(names: List[str] | None = None) -> Dict[str, Dict[str, Dict[int, float]]]:
    """Load CAID3 predictions, with metapredict taken from the 3.1.0 predictions.

    Parameters
    ----------
    names : list of str or None, optional
        Predictor names to load (CAID3 file stems). ``None`` loads every predictor.

    Returns
    -------
    dict
        Predictor name -> ``{target: {residue index: score}}``.
    """
    paths = dict(list_predictor_files(str(shared_data_path("caid3"))))
    for name in METAPREDICT_NAMES:
        paths[name] = str(metapredict_caid_path("caid3", name[-2:].upper()))
    if names is not None:
        missing = [n for n in names if n not in paths]
        assert not missing, f"no CAID3 predictions for {missing}"
        paths = {n: paths[n] for n in names}
    return {name: parse_caid(path) for name, path in paths.items()}


def rank_predictors(reference, predictions, dataset: str = DATASET_NAME) -> List[PredictorResult]:
    """Score every disorder predictor and return results sorted by AUC (best first).

    Parameters
    ----------
    reference : dict
        Target ID -> ``(sequence, labels)`` with labels 1 / 0 / -1 (excluded).
    predictions : dict
        Output of :func:`load_predictions`.
    dataset : str, optional
        Name recorded on each result (default ``"Disorder-PDB"``).

    Returns
    -------
    list of PredictorResult
        Disorder predictors only (binding/linker flavours excluded), best AUC first.
    """
    ranking = []
    for name, preds in predictions.items():
        result = evaluate_predictor(name, dataset, reference, preds)
        if result is not None and is_disorder_predictor(name):
            ranking.append(result)
    ranking.sort(key=lambda r: r.auc, reverse=True)
    return ranking


def print_metapredict_ranks(results: List[PredictorResult]) -> None:
    """Print each metapredict version's value and rank on every metric."""
    for attr in ["auc", "aps", "fmax", "mcc"]:
        order = sorted(results, key=lambda r: getattr(r, attr), reverse=True)
        names = [r.name for r in order]
        summary = ", ".join(
            f"{label.split()[-1]} {getattr(order[names.index(name)], attr):.3f} "
            f"(rank {names.index(name) + 1}/{len(names)})"
            for name, label in METAPREDICT_NAMES.items())
        print(f"{attr.upper():5s} {summary}")


# --------------------------------------------------------------------------- #
# Colours
# --------------------------------------------------------------------------- #
def predictor_colors(names: List[str]) -> Dict[str, str]:
    """Fixed colour for each predictor: metapredict, then the top CAID3 methods."""
    colors, extra = {}, 0
    for name in names:
        if name in METAPREDICT_COLORS:
            colors[name] = METAPREDICT_COLORS[name]
        elif name in CAID3_PREDICTOR_COLORS:
            colors[name] = CAID3_PREDICTOR_COLORS[name]
        else:
            colors[name] = EXTRA_COLORS[extra % len(EXTRA_COLORS)]
            extra += 1
    return colors


def _line_style(name: str) -> Tuple[float, int]:
    """(linewidth, zorder): emphasise metapredict, V3 most of all (widths as in Figure S11)."""
    if name == "Metapredict-v3":
        return 2.0, 7
    if name in METAPREDICT_NAMES:
        return 1.3, 6
    return 1.0, 4


# --------------------------------------------------------------------------- #
# Ranking bars
# --------------------------------------------------------------------------- #
def metric_panel(ax, ranking: List[PredictorResult], attr: str, ylabel: str,
                 label: str, top_n: int = 40) -> None:
    """Draw a CAID-style ranking of the top ``top_n`` predictors by one metric.

    Bars are black, ordered high to low; metapredict versions are recoloured. The
    y-axis is zoomed to the plotted range so differences among the top methods are
    visible. ``label`` and ``n`` (bars shown) are written inside the axes.
    """
    order = sorted(ranking, key=lambda r: getattr(r, attr), reverse=True)[:top_n]
    names = [display_name(r.name) for r in order]
    scores = [getattr(r, attr) for r in order]
    x = np.arange(1, len(scores) + 1)

    ax.bar(x, scores, color=FIELD_COLOR)
    for r, xi, score in zip(order, x, scores):
        if r.name in METAPREDICT_COLORS:
            ax.bar(xi, score, color=METAPREDICT_COLORS[r.name], linewidth=0.5,
                   edgecolor="black")

    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xticks(x, names, fontsize=7, rotation=45, horizontalalignment="right",
                  rotation_mode="anchor")
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlim(0.5, len(scores) + 0.5)
    lo, hi = min(scores), max(scores)
    span = (hi - lo) or 1.0
    ax.set_ylim(max(0.0, lo - 0.05 * span), hi + 0.15 * span)
    ax.text(0.99, 0.95, f"{label}  ($n = {len(names)}$)", transform=ax.transAxes,
            ha="right", va="top", fontsize=8)


# --------------------------------------------------------------------------- #
# ROC / PR against the field
# --------------------------------------------------------------------------- #
def highlighted(ranking: List[PredictorResult]) -> List[PredictorResult]:
    """Top-N predictors by AUC plus every metapredict version."""
    top = list(ranking[:N_TOP])
    present = {r.name for r in top}
    for r in ranking:
        if r.name in METAPREDICT_NAMES and r.name not in present:
            top.append(r)
            present.add(r.name)
    return top


def positive_rate(reference) -> float:
    """Fraction of scored residues that are disordered (the PR-curve baseline)."""
    pos = tot = 0
    for _seq, lab in reference.values():
        scored = lab[lab >= 0]
        pos += int((scored == 1).sum())
        tot += int(scored.size)
    return pos / tot


def field_curve_panel(ax, kind: str, reference, predictions,
                      ranking: List[PredictorResult], legend: bool = True) -> None:
    """ROC (``kind='roc'``) or precision-recall (``kind='pr'``) curves against the field.

    The top-N predictors (by AUC) and metapredict V1/V2/V3 are drawn in their fixed
    colours; every other disorder predictor is drawn unlabelled in light grey. The
    legend names predictors only (no metric values, to keep text off the panel);
    each highlighted predictor's AUC or APS and rank are printed instead.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    kind : str
        ``"roc"`` (true against false positive rate) or ``"pr"`` (precision against
        recall).
    reference : dict
        Target ID -> ``(sequence, labels)`` with labels 1 / 0 / -1 (excluded).
    predictions : dict
        Output of :func:`load_predictions`.
    ranking : list of PredictorResult
        Output of :func:`rank_predictors` (best AUC first).
    legend : bool, optional
        Draw the predictor legend in the panel's empty corner. Figure 2 draws it in B
        only (below the PR baseline, where no curve runs), since A uses the same
        colours. The PR legend is checked to sit below the baseline using the layout
        at call time, so give the figure ``layout="constrained"`` (or its final
        layout) rather than calling ``tight_layout`` afterwards.
    """
    assert kind in ("roc", "pr")
    attr = "auc" if kind == "roc" else "aps"
    curve = roc_curve_xy if kind == "roc" else pr_curve_xy
    series = highlighted(ranking)
    series_names = {r.name for r in series}

    # background: every other predictor, unlabelled
    n_field = 0
    for r in ranking:
        if r.name in series_names:
            continue
        y_true, y_score, _, _ = collect_pooled(reference, predictions[r.name])
        ax.plot(*curve(y_true, y_score), color=FIELD_CURVE_COLOR, lw=0.6, zorder=1)
        n_field += 1
    print(f"{n_field} other predictors drawn in grey")

    # foreground: top predictors and metapredict
    by_metric = sorted(ranking, key=lambda r: getattr(r, attr), reverse=True)
    rank_of = {r.name: i + 1 for i, r in enumerate(by_metric)}
    colors = predictor_colors([r.name for r in series])
    for r in series:
        y_true, y_score, _, _ = collect_pooled(reference, predictions[r.name])
        lw, z = _line_style(r.name)
        ax.plot(*curve(y_true, y_score), color=colors[r.name], lw=lw, zorder=z,
                label=display_name(r.name))
        print(f"{display_name(r.name):16s} {attr.upper()} {getattr(r, attr):.3f} "
              f"(rank {rank_of[r.name]}/{len(ranking)})")
    ax.plot([], [], color=FIELD_CURVE_COLOR, lw=1.0, label=FIELD_LABEL)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    legend_style = {"fontsize": LEGEND_FONTSIZE, "handlelength": 1.4, "labelspacing": 0.3}
    if kind == "roc":
        ax.plot([0, 1], [0, 1], ls=":", color="black", lw=0.6, zorder=2)
        ax.set_xlabel("False positive rate", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel("True positive rate", fontsize=LABEL_FONTSIZE)
        if legend:
            # the ROC curves leave the lower right of the panel mostly empty
            ax.legend(loc="lower right", **legend_style)
    else:
        baseline = positive_rate(reference)
        ax.axhline(baseline, ls=":", color="black", lw=0.6, zorder=2)
        ax.set_xlabel("Recall", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel("Precision", fontsize=LABEL_FONTSIZE)
        if legend:
            # no PR curve runs below the random baseline, so the legend goes there in
            # two columns (top CAID3 methods, then metapredict and the field)
            leg = ax.legend(loc="lower left", ncol=2, columnspacing=1.0, borderaxespad=0.3,
                            **legend_style)
            ax.figure.draw_without_rendering()
            legend_top = ax.transData.inverted().transform(leg.get_window_extent().corners()[1])[1]
            assert legend_top < baseline, (
                f"PR legend reaches precision {legend_top:.2f}, above the baseline {baseline:.2f}")


# --------------------------------------------------------------------------- #
# Per-residue traces
# --------------------------------------------------------------------------- #
def _runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return ``(first, last)`` 0-based indices of every contiguous run of True in ``mask``."""
    if not mask.any():
        return []
    idx = np.flatnonzero(mask)
    breaks = np.flatnonzero(np.diff(idx) != 1)
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    return list(zip(starts.tolist(), ends.tolist()))


def label_track(ax, labels: np.ndarray) -> None:
    """Draw one protein's Disorder-PDB labels as a thin bar (as in Figures 2D/E and S10).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        A short axes sharing the residue axis of the trace above it.
    labels : numpy.ndarray
        Per-residue reference labels (1 disordered, 0 ordered, -1 unscored); index
        ``i`` is residue ``i + 1`` and spans x from ``i + 0.5`` to ``i + 1.5``.
    """
    for state, color in DISORDER_PDB_COLORS.items():
        hatch = "///" if state == -1 else None
        for first, last in _runs(labels == state):
            ax.axvspan(first + 0.5, last + 1.5, facecolor=color, lw=0, hatch=hatch,
                       edgecolor=HATCH_EDGE_COLOR if hatch else "none")
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def label_track_handles() -> List[Patch]:
    """Legend entries for the Disorder-PDB bar (disordered, ordered, unscored)."""
    return [Patch(facecolor=color, hatch="///" if state == -1 else None,
                  edgecolor=HATCH_EDGE_COLOR if state == -1 else "none",
                  label=DISORDER_PDB_LABELS[state])
            for state, color in DISORDER_PDB_COLORS.items()]


def trace_panel_height(n_targets: int) -> float:
    """Height in inches of a :func:`trace_panel` figure with ``n_targets`` proteins."""
    block = TRACE_HEIGHT_IN + TRACE_TO_BAR_GAP_IN + LABEL_BAR_HEIGHT_IN
    return n_targets * block + (n_targets - 1) * BLOCK_GAP_IN


def trace_panel(fig, targets: List[str], series: List[str], reference, predictions,
                legend_predictors: bool = True) -> List:
    """Per-residue disorder traces over several proteins, each above its Disorder-PDB bar.

    Each protein gets a trace axes (every predictor in ``series`` in its fixed colour,
    with the 0.5 threshold dashed) and, directly below it, a thin bar of the
    Disorder-PDB labels (orange disordered, green ordered, hatched grey unscored). The
    protein's DisProt ID sits above its trace; one "Disorder score" label serves every
    trace; "Residue" labels the bottom axis only.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to draw into, ``trace_panel_height(len(targets))`` inches tall (the
        tick labels, axis labels and legend sit outside it and are kept by the
        tight bounding box used when saving).
    targets : list of str
        DisProt IDs of the proteins, top to bottom.
    series : list of str
        Predictor names (CAID3 file stems), drawn in their fixed colours.
    reference : dict
        Target ID -> ``(sequence, labels)``.
    predictions : dict
        Output of :func:`load_predictions`.
    legend_predictors : bool, optional
        Include the predictor lines in the legend. Figure 2C leaves them out because
        Figure 2B keys the same colours; the Disorder-PDB key is always drawn.

    Returns
    -------
    list of matplotlib.axes.Axes
        The trace axes, one per target, top first.
    """
    colors = predictor_colors(series)
    n_targets = len(targets)
    height_in = trace_panel_height(n_targets)
    # rows per protein: trace, gap, bar, then a gap before the next protein
    ratios: List[float] = []
    for i in range(n_targets):
        ratios += [TRACE_HEIGHT_IN, TRACE_TO_BAR_GAP_IN, LABEL_BAR_HEIGHT_IN]
        if i < n_targets - 1:
            ratios.append(BLOCK_GAP_IN)
    grid = GridSpec(len(ratios), 1, figure=fig, height_ratios=ratios, hspace=0,
                    left=0.0, right=1.0, top=1.0, bottom=0.0)

    trace_axes = []
    for i, tid in enumerate(targets):
        seq, labels = reference[tid]
        ax = fig.add_subplot(grid[4 * i, 0])
        ax_bar = fig.add_subplot(grid[4 * i + 2, 0], sharex=ax)
        for name in series:
            preds = predictions[name].get(tid)
            assert preds, f"{name} has no prediction for {tid}"
            lw, z = _line_style(name)
            ax.plot(*target_trace(labels, preds), color=colors[name], lw=lw, zorder=z)
        ax.axhline(DISORDER_THRESHOLD, ls="--", color="black", lw=0.5, zorder=2)
        ax.set_xlim(1, len(seq))
        ax.set_ylim(-0.02, 1.02)
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(labelsize=TICK_FONTSIZE, labelbottom=False)
        ax.text(0.0, 1.03, tid, transform=ax.transAxes, fontsize=TICK_FONTSIZE,
                ha="left", va="bottom")
        label_track(ax_bar, labels)
        ax_bar.tick_params(labelsize=TICK_FONTSIZE)
        trace_axes.append(ax)
    ax_bar.set_xlabel("Residue", fontsize=LABEL_FONTSIZE)
    # one y label for the whole stack, left of the tick labels
    fig.supylabel("Disorder score", x=-YLABEL_OFFSET_IN / fig.get_figwidth(), ha="right",
                  fontsize=LABEL_FONTSIZE)

    legend_style = {"fontsize": LEGEND_FONTSIZE, "handlelength": 1.4, "columnspacing": 1.2,
                    "borderaxespad": 0, "borderpad": 0}
    handles = label_track_handles()
    if not legend_predictors:
        # a short legend shares the top line with the first protein's ID, right-aligned
        header_base = 0.03 * TRACE_HEIGHT_IN / height_in
        fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.0 + header_base),
                   ncol=len(handles), **legend_style)
        return trace_axes
    # a long legend goes above the first protein's ID
    handles = [plt.Line2D([0], [0], color=colors[n], lw=_line_style(n)[0],
                          label=display_name(n)) for n in series] + handles
    header_top = (0.03 * TRACE_HEIGHT_IN + HEADER_HEIGHT_IN) / height_in
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0 + header_top),
               ncol=LEGEND_COLUMNS, **legend_style)
    return trace_axes
