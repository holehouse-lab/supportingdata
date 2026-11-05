"""Core library for reproducing the CAID3 disorder evaluation locally.

This module parses the CAID Disorder-PDB reference FASTA file
(``data/disorder_pdb.fasta``) and the per-predictor ``.caid`` prediction files,
then computes the same residue-level metrics reported in CAID3: the area under the
ROC curve (AUC), the average precision score (APS), the maximum F1 score
(``Fmax``), the Matthews correlation coefficient at the ``Fmax`` threshold, and
the fraction of reference targets a method covers.

The scoring conventions follow the CAID assessment:

* Reference labels use ``1`` for disordered (positive), ``0`` for ordered
  (negative), and ``-`` for residues excluded from scoring.
* A predictor's continuous propensity score (third column of a ``.caid`` file)
  is used when present. When that field is empty, the binary state (final
  column) is used as a proxy, matching the CAID rule that "if scores are
  unavailable, states are used as proxies".
* Predictions are aligned to the reference by residue index rather than by
  positional order, so missing residues or off-by-one records are handled
  gracefully.

Metrics are computed on the pooled set of covered, non-excluded residues, which
reproduces the dataset-level numbers reported in the CAID3 Disorder-PDB ranking
table.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
)


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def parse_reference(path: str) -> Dict[str, Tuple[str, np.ndarray]]:
    """Parse a CAID reference FASTA file.

    Each entry spans three lines: a header (``>ID``), the amino-acid sequence,
    and a label string of the same length using ``1`` (disordered), ``0``
    (ordered), and ``-`` (excluded).

    Parameters
    ----------
    path : str
        Path to the reference FASTA file.

    Returns
    -------
    dict
        Mapping of target ID to a tuple ``(sequence, labels)`` where ``labels``
        is an integer array with ``1``/``0`` for scored residues and ``-1`` for
        excluded residues.
    """
    entries: Dict[str, Tuple[str, np.ndarray]] = {}
    with open(path) as handle:
        lines = [ln.rstrip("\n") for ln in handle]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            target_id = lines[i][1:].strip()
            seq = lines[i + 1]
            label_str = lines[i + 2]
            labels = np.array(
                [-1 if c == "-" else int(c) for c in label_str], dtype=np.int8
            )
            entries[target_id] = (seq, labels)
            i += 3
        else:
            i += 1
    return entries


def parse_caid(path: str) -> Dict[str, Dict[int, float]]:
    """Parse a ``.caid`` prediction file into per-target index->score maps.

    Parameters
    ----------
    path : str
        Path to a ``.caid`` file. Records are ``index<TAB>aa<TAB>score`` with an
        optional trailing binary ``state`` column.

    Returns
    -------
    dict
        Mapping of target ID to a dict of ``{residue_index: score}``. When the
        score field is empty, the binary state column is used as the score.
    """
    predictions: Dict[str, Dict[int, float]] = {}
    current: Optional[str] = None
    rows: Dict[int, float] = {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    predictions[current] = rows
                current = line[1:].strip()
                rows = {}
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                idx = int(parts[0])
            except ValueError:
                continue
            score_field = parts[2].strip()
            if score_field == "":
                # Fall back to the binary state (last column) as a proxy score.
                state_field = parts[-1].strip()
                try:
                    score = float(state_field)
                except ValueError:
                    continue
            else:
                try:
                    score = float(score_field)
                except ValueError:
                    continue
            rows[idx] = score
        if current is not None:
            predictions[current] = rows
    return predictions


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
@dataclass
class PredictorResult:
    """Container for a predictor's pooled metrics on one reference set."""

    name: str
    dataset: str
    auc: float
    aps: float
    fmax: float
    mcc: float          # MCC at the Fmax threshold
    mcc_opt: float      # maximum MCC over all thresholds (the MCC analog of Fmax)
    fmax_threshold: float
    coverage: float
    n_targets_covered: int
    n_targets_total: int
    n_residues: int


def _fmax_and_mcc(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[float, float, float]:
    """Compute Fmax, the MCC at the Fmax threshold, and that threshold.

    ``Fmax`` is the maximum F1 score over all decision thresholds along the
    precision-recall curve. The Matthews correlation coefficient is reported at
    the threshold that maximises F1.
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    # precision/recall have length len(thresholds) + 1; drop the trailing point
    # that corresponds to recall == 0 and has no associated threshold.
    p = precision[:-1]
    r = recall[:-1]
    denom = p + r
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.where(denom > 0, 2 * p * r / denom, 0.0)
    if f1.size == 0:
        return 0.0, 0.0, 0.5
    best = int(np.argmax(f1))
    fmax = float(f1[best])
    thr = float(thresholds[best])

    pred = (y_score >= thr).astype(np.int8)
    mcc = _mcc(y_true, pred)
    return fmax, mcc, thr


def _mcc_opt(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Maximum Matthews correlation coefficient over all decision thresholds.

    Computed by sweeping every distinct score cut-off (the confusion matrix is
    accumulated from the score-sorted labels), so it is the threshold-optimised
    MCC — the direct analog of Fmax for F1.
    """
    if y_true.size == 0 or len(np.unique(y_true)) < 2:
        return 0.0
    order = np.argsort(y_score, kind="mergesort")[::-1]
    y = y_true[order].astype(np.int64)
    s = y_score[order]
    # Threshold at the last index of each run of equal scores.
    idx = np.r_[np.flatnonzero(np.diff(s)), y.size - 1]
    tp = np.cumsum(y)[idx].astype(np.float64)
    fp = (1 + idx - np.cumsum(y)[idx]).astype(np.float64)
    pos, neg = tp[-1], fp[-1]
    fn = pos - tp
    tn = neg - fp
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    with np.errstate(divide="ignore", invalid="ignore"):
        mcc = np.where(denom > 0, (tp * tn - fp * fn) / denom, 0.0)
    return float(np.max(mcc)) if mcc.size else 0.0


def _mcc(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Matthews correlation coefficient from binary vectors."""
    tp = float(np.sum((y_pred == 1) & (y_true == 1)))
    tn = float(np.sum((y_pred == 0) & (y_true == 0)))
    fp = float(np.sum((y_pred == 1) & (y_true == 0)))
    fn = float(np.sum((y_pred == 0) & (y_true == 1)))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def collect_pooled(
    reference: Dict[str, Tuple[str, np.ndarray]],
    predictions: Dict[str, Dict[int, float]],
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Pool covered, non-excluded residues across all reference targets.

    Parameters
    ----------
    reference : dict
        Output of :func:`parse_reference` for one dataset.
    predictions : dict
        Output of :func:`parse_caid` for one predictor.

    Returns
    -------
    tuple
        ``(y_true, y_score, n_targets_covered, n_targets_total)`` where the
        arrays contain one entry per scored residue that the predictor covered.
    """
    y_true: List[int] = []
    y_score: List[float] = []
    covered = 0
    for target_id, (_seq, labels) in reference.items():
        pred = predictions.get(target_id)
        if not pred:
            continue
        covered += 1
        for pos, label in enumerate(labels, start=1):  # reference is 1-indexed
            if label < 0:
                continue
            score = pred.get(pos)
            if score is None:
                continue
            y_true.append(int(label))
            y_score.append(score)
    return (
        np.asarray(y_true, dtype=np.int8),
        np.asarray(y_score, dtype=np.float64),
        covered,
        len(reference),
    )


def evaluate_predictor(
    name: str,
    dataset: str,
    reference: Dict[str, Tuple[str, np.ndarray]],
    predictions: Dict[str, Dict[int, float]],
) -> Optional[PredictorResult]:
    """Evaluate one predictor against one reference set.

    Returns ``None`` when the predictor covers no residues or the reference has
    only a single class among the covered residues (metrics undefined).
    """
    y_true, y_score, covered, total = collect_pooled(reference, predictions)
    if y_true.size == 0 or len(np.unique(y_true)) < 2:
        return None
    auc = float(roc_auc_score(y_true, y_score))
    aps = float(average_precision_score(y_true, y_score))
    fmax, mcc, thr = _fmax_and_mcc(y_true, y_score)
    mcc_opt = _mcc_opt(y_true, y_score)
    return PredictorResult(
        name=name,
        dataset=dataset,
        auc=auc,
        aps=aps,
        fmax=fmax,
        mcc=mcc,
        mcc_opt=mcc_opt,
        fmax_threshold=thr,
        coverage=covered / total if total else 0.0,
        n_targets_covered=covered,
        n_targets_total=total,
        n_residues=int(y_true.size),
    )


def list_predictor_files(pred_dir: str) -> List[Tuple[str, str]]:
    """Return ``(name, path)`` for every ``.caid`` file in ``pred_dir``."""
    files = sorted(glob.glob(os.path.join(pred_dir, "*.caid")))
    return [(os.path.basename(f)[:-5], f) for f in files]


def roc_curve_xy(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(fpr, tpr)`` for a ROC curve."""
    from sklearn.metrics import roc_curve

    fpr, tpr, _ = roc_curve(y_true, y_score)
    return fpr, tpr


def pr_curve_xy(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(recall, precision)`` for a precision-recall curve."""
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    return recall, precision


def target_trace(
    labels: np.ndarray, predictions: Dict[int, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-position scores aligned to a reference target.

    Parameters
    ----------
    labels : np.ndarray
        Reference label array for the target (from :func:`parse_reference`).
    predictions : dict
        ``{residue_index: score}`` map for the target (may be missing residues).

    Returns
    -------
    tuple
        ``(positions, scores)`` where ``positions`` is 1-based and ``scores``
        holds the predicted propensity (``np.nan`` where the method has no value).
    """
    positions = np.arange(1, len(labels) + 1)
    scores = np.array(
        [predictions.get(int(p), np.nan) for p in positions], dtype=np.float64
    )
    return positions, scores
