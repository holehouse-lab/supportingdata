"""CAID3 Disorder-PDB scored separately for disordered regions of different lengths.

For each IDR-length bin, the reference keeps every ordered residue but only the
disordered residues that lie in an IDR (a maximal run of consecutive residues
labelled disordered) whose length falls in the bin; every other disordered residue
is excluded. Each bin is then scored exactly as Figure S3 scores the full set
(``caid3_panels.rank_predictors``, i.e. the unchanged CAID3 scoring code), so every
method's Fmax threshold is chosen within each bin.

The paired bootstrap re-implements the pooled metrics with per-target weights,
so that 1,000 resamples take seconds rather than an hour; ``check_fast_metrics``
confirms that it reproduces the CAID3 scoring code exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

from caid3_panels import (
    METAPREDICT_COLORS,
    METAPREDICT_NAMES,
    FIELD_CURVE_COLOR,
    display_name,
    metric_panel,
    rank_predictors,
)
from caid_eval import PredictorResult
from figtools import CAID3_PREDICTOR_COLORS, V3_COLOR

# IDR-length bins: (shortest IDR, longest IDR or None for no limit, axis label), in residues.
IDR_LENGTH_BINS: List[Tuple[int, Optional[int], str]] = [
    (1, 19, "<20"),
    (20, 49, "20–49"),
    (50, 99, "50–99"),
    (100, None, "≥100"),
]

METRIC_LABELS = {"auc": "AUC", "aps": "APS", "fmax": "Fmax", "mcc": "MCC"}
REFERENCE_PREDICTOR = "Metapredict-v3"  # every bootstrap comparison is against V3

# Paired bootstrap over targets, as in the CAID assessments. The seed is arbitrary
# but fixed so the reported intervals are reproducible.
N_BOOTSTRAP = 1000
BOOTSTRAP_SEED = 20261030
CONFIDENCE_PERCENT = 95.0

# Tolerance for the fast metrics against the CAID3 scoring code (float round-off only).
FAST_METRIC_TOLERANCE = 1e-9

# The highlighted CAID3 methods are the top five on the full Disorder-PDB set
# (by AUC), drawn in the same fixed colours as Figure 2.
HIGHLIGHTED = list(CAID3_PREDICTOR_COLORS)


# --------------------------------------------------------------------------- #
# Stratifying the reference
# --------------------------------------------------------------------------- #
def idr_run_lengths(labels: np.ndarray) -> np.ndarray:
    """Length of the IDR that each residue belongs to.

    Parameters
    ----------
    labels : numpy.ndarray
        Per-residue reference labels: 1 disordered, 0 ordered, -1 excluded.

    Returns
    -------
    numpy.ndarray
        For residues labelled disordered, the length (residues) of the maximal run of
        consecutive disordered residues containing them; 0 for every other residue.
    """
    lengths = np.zeros(labels.size, dtype=np.int64)
    padded = np.r_[0, (labels == 1).astype(np.int8), 0]
    edges = np.flatnonzero(np.diff(padded))
    for start, end in zip(edges[::2], edges[1::2]):
        lengths[start:end] = end - start
    return lengths


def in_bin(lengths: np.ndarray, shortest: int, longest: Optional[int]) -> np.ndarray:
    """Boolean mask of IDR lengths inside ``[shortest, longest]`` (``longest=None``: no limit)."""
    upper = np.inf if longest is None else longest
    return (lengths >= shortest) & (lengths <= upper)


def stratify_reference(reference: Dict[str, Tuple[str, np.ndarray]], shortest: int,
                       longest: Optional[int]) -> Dict[str, Tuple[str, np.ndarray]]:
    """Reference with disordered residues restricted to IDRs of one length bin.

    Parameters
    ----------
    reference : dict
        Target ID -> ``(sequence, labels)``, labels 1 / 0 / -1.
    shortest, longest : int, int or None
        IDR length bin in residues (inclusive); ``longest=None`` means no upper limit.

    Returns
    -------
    dict
        Same targets; ordered residues unchanged, disordered residues outside the bin
        relabelled -1 (excluded).
    """
    stratified = {}
    for target, (sequence, labels) in reference.items():
        new_labels = labels.copy()
        outside = (labels == 1) & ~in_bin(idr_run_lengths(labels), shortest, longest)
        new_labels[outside] = -1
        stratified[target] = (sequence, new_labels)
    return stratified


def describe_bins(reference: Dict[str, Tuple[str, np.ndarray]]) -> None:
    """Print the number of IDRs and disordered residues in each bin (ordered residues are shared)."""
    n_ordered = sum(int((labels == 0).sum()) for _seq, labels in reference.values())
    print(f"every bin keeps all {n_ordered} ordered residues")
    for shortest, longest, label in IDR_LENGTH_BINS:
        n_idrs = n_residues = 0
        for _seq, labels in reference.values():
            lengths = idr_run_lengths(labels)
            keep = (labels == 1) & in_bin(lengths, shortest, longest)
            n_residues += int(keep.sum())
            starts = np.flatnonzero(np.diff(np.r_[0, keep.astype(np.int8)]) == 1)
            n_idrs += starts.size
        print(f"IDRs {label:>6s} aa: {n_idrs:4d} IDRs, {n_residues:6d} disordered residues")


def rank_by_idr_length(reference: Dict[str, Tuple[str, np.ndarray]],
                       predictions: dict) -> Dict[str, List[PredictorResult]]:
    """Score every CAID3 disorder predictor on each IDR-length bin.

    Parameters
    ----------
    reference : dict
        The full CAID3 Disorder-PDB reference.
    predictions : dict
        Output of ``caid3_panels.load_predictions``.

    Returns
    -------
    dict
        Bin label -> list of ``PredictorResult`` (best AUC first), as for Figure S3.
    """
    rankings = {}
    for shortest, longest, label in IDR_LENGTH_BINS:
        stratified = stratify_reference(reference, shortest, longest)
        rankings[label] = rank_predictors(stratified, predictions, dataset=f"IDRs {label} aa")
    return rankings


def rank_of(ranking: List[PredictorResult], name: str, attr: str) -> int:
    """1-based rank of predictor ``name`` by metric ``attr`` (higher is better)."""
    order = sorted(ranking, key=lambda r: getattr(r, attr), reverse=True)
    return [r.name for r in order].index(name) + 1


def print_metapredict_ranks_by_bin(rankings: Dict[str, List[PredictorResult]]) -> None:
    """Print each metapredict version's value and rank on every metric in every bin."""
    for label, ranking in rankings.items():
        by_name = {r.name: r for r in ranking}
        print(f"IDRs {label} aa ({len(ranking)} predictors)")
        for attr, metric in METRIC_LABELS.items():
            cells = ", ".join(f"{version.split()[-1]} {getattr(by_name[name], attr):.3f} "
                              f"(rank {rank_of(ranking, name, attr)})"
                              for name, version in METAPREDICT_NAMES.items())
            print(f"  {metric:5s} {cells}")


# --------------------------------------------------------------------------- #
# Fast pooled metrics and the paired bootstrap
# --------------------------------------------------------------------------- #
@dataclass
class SortedScores:
    """One predictor's scored residues, sorted once by descending score.

    ``target_index`` gives each residue's target (an index into the target list) so a
    bootstrap resample is just a per-target weight; ``tie_ends`` marks the last residue
    of each run of equal scores (one decision threshold per run).
    """

    labels: np.ndarray
    target_index: np.ndarray
    tie_ends: np.ndarray


def sort_scores(reference: Dict[str, Tuple[str, np.ndarray]], predictions: Dict[str, Dict[int, float]],
                targets: List[str]) -> SortedScores:
    """Collect the covered, non-excluded residues of one predictor and sort them by score.

    Uses the same residues as ``caid_eval.collect_pooled``.
    """
    labels, scores, target_index = [], [], []
    for i, target in enumerate(targets):
        prediction = predictions.get(target)
        if not prediction:
            continue
        for position, label in enumerate(reference[target][1], start=1):  # reference is 1-indexed
            if label < 0 or position not in prediction:
                continue
            labels.append(int(label))
            scores.append(prediction[position])
            target_index.append(i)
    score_array = np.asarray(scores, dtype=np.float64)
    order = np.argsort(-score_array, kind="mergesort")
    sorted_scores = score_array[order]
    tie_ends = np.r_[np.flatnonzero(np.diff(sorted_scores)), sorted_scores.size - 1]
    return SortedScores(np.asarray(labels, dtype=np.float64)[order],
                        np.asarray(target_index, dtype=np.int64)[order], tie_ends)


def weighted_metrics(scores: SortedScores, target_weights: np.ndarray) -> Dict[str, float]:
    """AUC, APS, Fmax and MCC at the Fmax threshold, with each target counted ``weight`` times.

    With all weights equal to 1 this is the CAID3 pooled scoring; with bootstrap
    multiplicities it equals scoring the concatenated resample.
    """
    weights = target_weights[scores.target_index].astype(np.float64)
    true_pos = np.cumsum(weights * scores.labels)[scores.tie_ends]
    false_pos = np.cumsum(weights * (1.0 - scores.labels))[scores.tie_ends]
    n_pos, n_neg = true_pos[-1], false_pos[-1]
    auc = np.trapezoid(np.r_[0.0, true_pos / n_pos], np.r_[0.0, false_pos / n_neg])
    called = true_pos + false_pos
    precision = np.divide(true_pos, called, out=np.ones_like(called), where=called > 0)
    recall = true_pos / n_pos
    aps = np.sum(np.diff(np.r_[0.0, recall]) * precision)
    denominator = precision + recall
    f1 = np.divide(2 * precision * recall, denominator, out=np.zeros_like(denominator),
                   where=denominator > 0)
    best = int(np.argmax(f1))
    tp, fp = true_pos[best], false_pos[best]
    fn, tn = n_pos - tp, n_neg - fp
    mcc_denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / mcc_denominator if mcc_denominator > 0 else 0.0
    return {"auc": float(auc), "aps": float(aps), "fmax": float(f1[best]), "mcc": float(mcc)}


def check_fast_metrics(ranking: List[PredictorResult], sorted_by_name: Dict[str, SortedScores],
                       n_targets: int) -> None:
    """Assert that the fast metrics reproduce the CAID3 scoring code for every predictor given."""
    by_name = {r.name: r for r in ranking}
    for name, scores in sorted_by_name.items():
        fast = weighted_metrics(scores, np.ones(n_targets))
        for attr in METRIC_LABELS:
            difference = abs(fast[attr] - getattr(by_name[name], attr))
            assert difference < FAST_METRIC_TOLERANCE, f"{name} {attr}: fast metric off by {difference}"


def paired_bootstrap(reference: Dict[str, Tuple[str, np.ndarray]], predictions: dict,
                     ranking: List[PredictorResult], n_resamples: int = N_BOOTSTRAP,
                     seed: int = BOOTSTRAP_SEED) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """Paired bootstrap of every predictor that beats metapredict V3 on some metric.

    Targets are resampled with replacement (the same resamples for every predictor).
    In each resample a predictor and V3 are scored on the same drawn targets that
    the predictor covers, so predictors that skipped targets are compared with V3 on
    the targets they predicted.

    Parameters
    ----------
    reference : dict
        Reference (already stratified to one IDR-length bin).
    predictions : dict
        Output of ``caid3_panels.load_predictions``.
    ranking : list of PredictorResult
        ``rank_predictors`` result for this reference.
    n_resamples : int, optional
        Number of bootstrap resamples (default 1,000).
    seed : int, optional
        Random seed for the resamples.

    Returns
    -------
    dict
        Metric -> {predictor: (lower, upper)} central ``CONFIDENCE_PERCENT`` interval of
        (predictor - V3), for every predictor ahead of V3 on that metric.
    """
    targets = sorted(reference)
    by_name = {r.name: r for r in ranking}
    v3 = by_name[REFERENCE_PREDICTOR]
    ahead = {attr: [r.name for r in ranking if getattr(r, attr) > getattr(v3, attr)]
             for attr in METRIC_LABELS}
    to_test = sorted({name for names in ahead.values() for name in names})

    sorted_by_name = {name: sort_scores(reference, predictions[name], targets)
                      for name in to_test + [REFERENCE_PREDICTOR]}
    check_fast_metrics(ranking, sorted_by_name, len(targets))

    rng = np.random.default_rng(seed)
    multiplicities = [np.bincount(rng.integers(0, len(targets), len(targets)), minlength=len(targets))
                      for _ in range(n_resamples)]
    tail = (100.0 - CONFIDENCE_PERCENT) / 2
    intervals: Dict[str, Dict[str, Tuple[float, float]]] = {attr: {} for attr in METRIC_LABELS}
    for name in to_test:
        covered = np.zeros(len(targets))
        covered[np.unique(sorted_by_name[name].target_index)] = 1.0
        differences = {attr: [] for attr in METRIC_LABELS}
        for weights in multiplicities:
            other = weighted_metrics(sorted_by_name[name], weights)
            mine = weighted_metrics(sorted_by_name[REFERENCE_PREDICTOR], weights * covered)
            for attr in METRIC_LABELS:
                differences[attr].append(other[attr] - mine[attr])
        for attr in METRIC_LABELS:
            if name in ahead[attr]:
                lower, upper = np.percentile(differences[attr], [tail, 100.0 - tail])
                intervals[attr][name] = (float(lower), float(upper))
    return intervals


def significantly_ahead(intervals: Dict[str, Tuple[float, float]]) -> List[str]:
    """Predictors whose whole interval of (predictor - V3) lies above zero."""
    return [name for name, (lower, _upper) in intervals.items() if lower > 0]


# --------------------------------------------------------------------------- #
# Panels
# --------------------------------------------------------------------------- #
def _series_style(name: str) -> dict:
    """Line style for a highlighted predictor (metapredict V3 heaviest)."""
    if name == REFERENCE_PREDICTOR:
        return {"color": METAPREDICT_COLORS[name], "lw": 2.0, "ms": 4.0, "zorder": 7}
    if name in METAPREDICT_COLORS:
        return {"color": METAPREDICT_COLORS[name], "lw": 1.3, "ms": 3.0, "zorder": 6}
    return {"color": CAID3_PREDICTOR_COLORS[name], "lw": 1.0, "ms": 2.5, "zorder": 5}


def field_panel(ax, rankings: Dict[str, List[PredictorResult]], attr: str,
                show_legend: bool = False) -> None:
    """One metric across IDR-length bins for every CAID3 disorder predictor.

    Every predictor is a short light-grey tick in each bin; the top five CAID3 methods
    and metapredict V1/V2/V3 are drawn as lines across bins. metapredict V3's rank in
    each bin is written above the axes in the V3 colour.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    rankings : dict
        Output of :func:`rank_by_idr_length`.
    attr : str
        ``"auc"``, ``"aps"``, ``"fmax"`` or ``"mcc"``.
    show_legend : bool, optional
        Draw the legend (to the right of the axes).
    """
    labels = list(rankings)
    x = np.arange(len(labels))
    for xi, label in zip(x, labels):
        values = [getattr(r, attr) for r in rankings[label]]
        ax.plot(np.full(len(values), xi), values, linestyle="none", marker="_", ms=9, mew=0.8,
                color=FIELD_CURVE_COLOR, zorder=2)

    for name in HIGHLIGHTED + list(METAPREDICT_NAMES):
        values = [getattr(next(r for r in rankings[label] if r.name == name), attr) for label in labels]
        ax.plot(x, values, marker="o", label=display_name(name), **_series_style(name))

    top_of_axes = blended_transform_factory(ax.transData, ax.transAxes)
    for xi, label in zip(x, labels):
        ax.text(xi, 1.02, f"#{rank_of(rankings[label], REFERENCE_PREDICTOR, attr)}",
                transform=top_of_axes, ha="center", va="bottom", fontsize=7, color=V3_COLOR)

    ax.set_xticks(x, labels, fontsize=7)
    ax.set_xlim(-0.4, len(labels) - 0.6)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlabel("IDR length (residues)", fontsize=8)
    ax.set_ylabel(METRIC_LABELS[attr], fontsize=8)
    if show_legend:
        field_handle = Line2D([], [], linestyle="none", marker="_", ms=9, mew=0.8,
                              color=FIELD_CURVE_COLOR, label="other CAID3 predictors")
        handles, _ = ax.get_legend_handles_labels()
        ax.legend(handles=handles + [field_handle], fontsize=6.5, loc="upper left",
                  bbox_to_anchor=(1.02, 1.0), handlelength=1.6)


def ranking_panel(ax, ranking: List[PredictorResult], attr: str, label: str,
                  significant: List[str], top_n: int = 40) -> None:
    """Figure S3-style ranking bars with an asterisk over predictors significantly ahead of V3.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    ranking : list of PredictorResult
        One bin's ranking.
    attr : str
        Metric to rank by.
    label : str
        In-axis label (the bin).
    significant : list of str
        Predictors significantly ahead of metapredict V3 (from :func:`significantly_ahead`).
    top_n : int, optional
        Number of bars.
    """
    metric_panel(ax, ranking, attr, METRIC_LABELS[attr], label, top_n=top_n)
    order = sorted(ranking, key=lambda r: getattr(r, attr), reverse=True)[:top_n]
    for xi, result in enumerate(order, start=1):
        if result.name in significant:
            ax.text(xi, getattr(result, attr), "*", ha="center", va="bottom", fontsize=9)
