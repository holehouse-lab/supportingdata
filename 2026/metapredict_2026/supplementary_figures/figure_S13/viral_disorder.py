"""Loading and strip-plot helpers for Figure S13B/C (viral proteome disorder).

Taken from the original ``supplemental_figure_6_B_C.ipynb``. Changes: the
grouping function takes the lineage and disorder tables as arguments instead of
reading notebook globals, the lineage plot returns its figure instead of saving
it, and the table readers are factored out of the notebook cells. The jitter
functions use Python's / NumPy's global random state; the notebooks seed both.
"""

from __future__ import annotations

import random

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_single_gene_viruses(path: str) -> list[str]:
    """Return proteome IDs whose proteome contains a single protein.

    Parameters
    ----------
    path : str
        ``viral_proteome_id_to_count.tsv`` (proteome ID, organism ID, protein count).

    Returns
    -------
    list of str
        Proteome IDs with exactly one protein.
    """
    viruses_by_count = {}
    with open(path) as fh:
        for sl in fh.read().split('\n')[1:]:
            splitvals = sl.split('\t')
            if splitvals[0] != '':
                viruses_by_count[splitvals[0]] = int(splitvals[-1])
    return [v for v, count in viruses_by_count.items() if count == 1]


def load_lineages(path: str) -> dict[str, str]:
    """Return ``{proteome ID: comma-separated taxonomic lineage}``.

    Parameters
    ----------
    path : str
        ``viruses_w_lineage.tsv`` (UniProt proteome table, lineage in the last column).
    """
    viruses = {}
    with open(path) as fh:
        for sl in fh.read().split('\n')[1:]:
            splitvals = sl.split('\t')
            if splitvals[0] != '':
                viruses[splitvals[0]] = splitvals[-1]
    return viruses


def load_viral_disorder(path: str, exclude: list[str]) -> dict[str, dict[str, float]]:
    """Return per-proteome disorder, skipping proteomes listed in ``exclude``.

    Parameters
    ----------
    path : str
        ``viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv``.
    exclude : list of str
        Proteome IDs to drop (e.g. single-gene proteomes).

    Returns
    -------
    dict
        ``{proteome ID: {'per_protein': AFPP, 'total_proteome': TFAP}}``.
    """
    exclude = set(exclude)
    vir_dis_dict = {}
    with open(path) as fh:
        for line in fh.read().split('\n')[1:]:
            if line == '':
                continue
            splitvals = line.split('\t')
            if splitvals[0] not in exclude:
                vir_dis_dict[splitvals[0]] = {'per_protein': float(splitvals[1]),
                                              'total_proteome': float(splitvals[2])}
    return vir_dis_dict


def read_in_tsv(path_to_tsv: str) -> dict[str, list[float]]:
    """Read a host-split viral disorder table (proteome ID, AFPP, TFAP).

    Returns
    -------
    dict
        ``{'per_protein': [...], 'total_proteome': [...]}``.
    """
    per_protein = []
    total_proteome = []
    with open(path_to_tsv) as fh:
        lines = fh.read().split('\n')
    for line in lines[1:]:
        if line != '' and len(line.split('\t')) == 3:
            vals = line.split('\t')
            per_protein.append(float(vals[1]))
            total_proteome.append(float(vals[2]))
    return {'per_protein': per_protein, 'total_proteome': total_proteome}


# --------------------------------------------------------------------------- #
# Grouping by lineage
# --------------------------------------------------------------------------- #
def break_down_into_groups(groups, lineages, disorder, min_per_group=5, vals='per_protein'):
    """Collect disorder values for each lineage group with at least ``min_per_group`` proteomes.

    Parameters
    ----------
    groups : list of str
        Lineage names to collect (matched case-insensitively against each lineage).
    lineages : dict
        ``{proteome ID: lineage string}`` from :func:`load_lineages`.
    disorder : dict
        Per-proteome disorder from :func:`load_viral_disorder`.
    min_per_group : int, optional
        Groups with fewer proteomes are dropped.
    vals : {'per_protein', 'total_proteome'}
        Disorder metric (AFPP or TFAP).

    Returns
    -------
    dict
        ``{group: [disorder values]}`` in the order of ``groups``.
    """
    if vals not in ['total_proteome', 'per_protein']:
        raise Exception('vals must be total_proteome or per_protein')

    group_dict = {p: [] for p in groups}
    for a in lineages:
        t = [k.lower() for k in lineages[a].split(', ')]
        for v in groups:
            if v.lower() in t and a in disorder:
                group_dict[v].append(disorder[a][vals])

    return {v: group_dict[v] for v in group_dict if len(group_dict[v]) >= min_per_group}


# --------------------------------------------------------------------------- #
# Jittered strip plots
# --------------------------------------------------------------------------- #
def generate_uniformly_spaced_points(min_val, max_val, num_points):
    """Return ``num_points`` uniformly spaced values between ``min_val`` and ``max_val``."""
    if num_points < 2:
        raise ValueError("num_points must be at least 2 to generate a range of points.")
    step = (max_val - min_val) / (num_points - 1)
    return [min_val + i * step for i in range(num_points)]


def jitter_based_on_density(y_coords, grouping_by=0.01, min_val=0, max_val=0.5, base_jitter=0.15):
    """Jitter x positions so that denser stretches of y are spread wider.

    Parameters
    ----------
    y_coords : list of float
        Values to plot (y).
    grouping_by : float, optional
        Bin width in y used to estimate local density.
    min_val, max_val : float, optional
        x range of the strip; points are centred between them.
    base_jitter : float, optional
        Extra half-width added for bins with more than four points.

    Returns
    -------
    list
        ``[sorted y values, x positions]``.
    """
    data = sorted(y_coords)
    groups = []
    center = min_val + ((max_val - min_val) / 2)
    total_vals = len(y_coords)
    for i in np.arange(min(y_coords) - grouping_by, max(y_coords) + (grouping_by * 3), grouping_by):
        i = round(float(i), 5)
        totpoints = len([a for a in y_coords if a >= i and a < round(i + grouping_by, 5)])
        if totpoints == 0:
            continue
        if totpoints == 1:
            groups.extend([center])
            continue
        plus_minus = (totpoints / total_vals) / 2
        if totpoints <= 4:
            curvals = generate_uniformly_spaced_points(center - plus_minus, center + plus_minus, totpoints)
        elif totpoints < 10:
            curvals = generate_uniformly_spaced_points(center - (plus_minus + (base_jitter / 2)),
                                                       center + (plus_minus + (base_jitter / 2)), totpoints)
        else:
            curvals = generate_uniformly_spaced_points(center - (plus_minus + base_jitter),
                                                       center + (plus_minus + base_jitter), totpoints)
        random.shuffle(curvals)
        groups.extend(curvals)
    return [data, groups]


def graph_lineages(group_disorder_dict, figsize=(4.5, 4.5), ylabel='', marker='.', s=3, alpha=0.2):
    """Strip plot of disorder per lineage group with a black bar at each group mean.

    Parameters
    ----------
    group_disorder_dict : dict
        ``{group: [disorder values]}``, plotted left to right in dict order.
    figsize : tuple, optional
        Figure size (inches).
    ylabel : str, optional
        y-axis label.
    marker, s, alpha : optional
        Scatter marker, size and opacity.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, total = plt.subplots(ncols=1, nrows=1, figsize=figsize, dpi=150)
    for n, p in enumerate(group_disorder_dict):
        data = jitter_based_on_density(group_disorder_dict[p], min_val=n, max_val=n)
        total.scatter(data[1], data[0], marker=marker, s=s, alpha=alpha)
        avg = np.mean(data[0])
        total.plot([n - 0.2, n + 0.2], [avg, avg], color='black', linewidth=1)
    total.set_xticks(list(range(len(group_disorder_dict))), list(group_disorder_dict.keys()), rotation=90)
    total.set_xlim([-0.5, len(group_disorder_dict) - 0.25])
    total.set_ylabel(ylabel)
    plt.tight_layout()
    return fig


def generate_density_based_scatter_data(y_values, center=0, base_width=1.0, scale_factor=2.0):
    """Spread x positions around ``center`` in proportion to the local density of y.

    Parameters
    ----------
    y_values : list or numpy.ndarray
        Values to plot (y).
    center : float, optional
        Centre of the x distribution.
    base_width : float, optional
        x spread (normal s.d.) in the sparsest regions.
    scale_factor : float, optional
        Additional x spread at the densest region (density normalised to [0, 1]).

    Returns
    -------
    tuple
        ``(x_values, y_values)``.
    """
    kde = gaussian_kde(y_values)
    density = kde(y_values)
    density_normalized = (density - density.min()) / (density.max() - density.min())
    x_spread = base_width + scale_factor * density_normalized
    x_values = np.random.normal(loc=center, scale=x_spread, size=len(y_values))
    return x_values, y_values
