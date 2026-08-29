"""Shared helpers for the PENGUIN figure notebooks.

Everything in this module is a thin layer on top of the public
:class:`soursop.sssampling.SamplingQuality` API. The notebooks in this
directory import from here so that the figure code is written once and every
notebook stays focused on a single figure from Lotthammer & Holehouse (2025),
*Disentangling folding from energetic traps in simulations of disordered
proteins*, J. Chem. Inf. Model., doi:10.1021/acs.jcim.4c02005.

Data location
-------------
The trajectories behind these figures (~6 GB) live in the ``data/`` directory
next to this file, one folder per system with ``full/`` (full-Hamiltonian) and
``ev/`` (excluded-volume) replicates - see ``data/README.md``. They are too
large to ship in the SOURSOP repository; set the environment variable
``PENGUIN_DATA_DIR`` to point elsewhere if you keep them in another location.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms

from soursop.sssampling import SamplingQuality
from soursop.sstools import find_trajectory_files

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DEFAULT_DATA_DIR = Path(__file__).resolve().parent / "data"

# figures are written next to the notebooks
FIGURE_DIR = Path(__file__).resolve().parent / "figures"


def get_data_dir() -> Path:
    """Resolve the simulation data directory.

    Returns
    -------
    pathlib.Path
        ``$PENGUIN_DATA_DIR`` if set, otherwise :data:`DEFAULT_DATA_DIR`.

    Raises
    ------
    FileNotFoundError
        If the resolved directory does not exist, with a message explaining
        how to point the notebooks at the data.
    """
    data_dir = Path(os.environ.get("PENGUIN_DATA_DIR", DEFAULT_DATA_DIR)).expanduser()
    if not data_dir.is_dir():
        raise FileNotFoundError(
            f"PENGUIN simulation data not found at {data_dir}. Set the "
            "PENGUIN_DATA_DIR environment variable (or edit DEFAULT_DATA_DIR in "
            "penguin_helpers.py) to the directory containing the per-system "
            "folders described in data/README.md (fs_peptide/, polyGV_PGSK/, ...)."
        )
    return data_dir


def figure_path(*parts: str) -> Path:
    """Build (and create) an output path under the notebook ``figures/`` folder.

    Parameters
    ----------
    *parts : str
        Sub-directory and file name components, e.g. ``('fig4', 'fs.pdf')``.

    Returns
    -------
    pathlib.Path
        Absolute path; parent directories are created on demand.
    """
    out = FIGURE_DIR.joinpath(*parts)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


# ---------------------------------------------------------------------------
# Matplotlib defaults matching the manuscript
# ---------------------------------------------------------------------------


def set_manuscript_style() -> None:
    """Apply the rcParams used for the manuscript figures.

    Fonts are embedded as TrueType (type 42) so the PDFs stay editable in
    Illustrator / Affinity, axes lines are thin, and Arial is requested with
    a fallback to the default sans-serif family if it is not installed.
    """
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["axes.linewidth"] = 0.5
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


# ---------------------------------------------------------------------------
# Data discovery
# ---------------------------------------------------------------------------


def load_system(
    full_subdir: str,
    ev_subdir: str,
    num_replicates: int,
    method: str = "2D angle distributions",
    bwidth_degrees: float = 15.0,
    truncate: bool = False,
    full_kwargs: Optional[dict] = None,
    ev_kwargs: Optional[dict] = None,
    seed: int = 0,
    verbose: bool = True,
    **kwargs,
) -> SamplingQuality:
    """Discover replicate trajectories for a system and build a SamplingQuality.

    Parameters
    ----------
    full_subdir, ev_subdir : str
        Sub-directories of the data directory holding the full-Hamiltonian
        (FH) and excluded-volume (EV) replicate simulations respectively.
        Each is searched recursively for ``__traj.xtc`` / ``__START.pdb``
        pairs inside numbered replicate folders.
    num_replicates : int
        Upper bound on the replicate index accepted by
        :func:`soursop.sstools.find_trajectory_files`.
    method : str, optional
        ``'2D angle distributions'`` (joint phi/psi, the manuscript default)
        or ``'1D angle distributions'``.
    bwidth_degrees : float, optional
        Histogram bin width in degrees. The manuscript uses 15.
    truncate : bool, optional
        Truncate every trajectory to the shortest one. Needed whenever the FH
        and EV runs have different lengths.
    full_kwargs, ev_kwargs : dict, optional
        Extra keyword arguments for the two ``find_trajectory_files`` calls
        (e.g. ``exclude_dirs``).
    seed : int, optional
        Seed for the (unused here, but deterministic) EV reference sampler.
    verbose : bool, optional
        Print the truncation summary.
    **kwargs
        Forwarded to :class:`SamplingQuality`.

    Returns
    -------
    SamplingQuality
    """
    data_dir = get_data_dir()
    trj, top = find_trajectory_files(
        data_dir / full_subdir, num_replicates, **(full_kwargs or {})
    )
    ev_trj, ev_top = find_trajectory_files(
        data_dir / ev_subdir, num_replicates, **(ev_kwargs or {})
    )
    if len(trj) == 0 or len(ev_trj) == 0:
        raise FileNotFoundError(
            f"No replicates found under {data_dir / full_subdir} / {data_dir / ev_subdir}"
        )
    if verbose:
        print(f"{len(trj)} FH replicates, {len(ev_trj)} EV replicates")
    return SamplingQuality(
        trj,
        ev_trj,
        top_file=top[0],
        ref_top=ev_top[0],
        method=method,
        bwidth=np.deg2rad(bwidth_degrees),
        truncate=truncate,
        seed=seed,
        verbose=verbose,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# The PENGUIN four-panel figure
# ---------------------------------------------------------------------------


def _panel_label(fig, ax, label: str, fontsize: float, position=(-0.0825, 1.10)):
    trans = transforms.ScaledTranslation(-20 / 72, 7 / 72, fig.dpi_scale_trans)
    ax.text(
        position[0],
        position[1],
        label,
        transform=ax.transAxes + trans,
        fontsize=fontsize,
        fontweight="bold",
        va="top",
        ha="right",
    )


def penguin_figure(
    qual: SamplingQuality,
    dihedral: str = "2D",
    increment: int = 5,
    highlight: Optional[Tuple[float, float]] = None,
    divider: Optional[float] = None,
    color: str = "#7689de",
    figsize: Tuple[float, float] = (7, 6),
    dpi: int = 300,
    panel_labels: Sequence[str] = ("B", "C", "D", "E"),
    fontsize: float = 10,
    label_position: Tuple[float, float] = (-0.0825, 1.10),
    title: Optional[str] = None,
    ev_title: str = "Comparison to the Excluded Volume Limit",
    save_path: Optional[os.PathLike] = None,
):
    """Draw the four-panel PENGUIN summary used throughout the manuscript.

    The panels are, in reading order:

    1. per-residue Hellinger distance between each FH replicate and the EV
       reference (red points, one per replicate) with the replicate mean
       (black);
    2. per-residue all-to-all Hellinger distance between FH replicates;
    3. the spread (max - min across replicates) of panel 1;
    4. per-residue fractional helicity of the FH replicates.

    Parameters
    ----------
    qual : SamplingQuality
        A constructed sampling-quality object.
    dihedral : {'2D', 'phi', 'psi'}, optional
        Which comparison to draw. ``'2D'`` needs ``method='2D angle
        distributions'``; ``'phi'``/``'psi'`` need the 1D method.
    increment : int, optional
        Residue tick spacing.
    highlight : tuple of (float, float), optional
        Shade the residue range ``[x1, x2]`` in every panel (used for the
        fixed region of FS-peptide and the poorly sampled Ash1 region).
    divider : float, optional
        Draw a dashed vertical line at this residue position in every panel
        (used for the GV/PGSK boundary).
    color : str, optional
        Colour for the highlight / divider.
    figsize, dpi : optional
        Figure size in inches and resolution.
    panel_labels : sequence of str, optional
        Letters for the four panels. The manuscript reserves ``A`` for the
        sequence cartoon drawn separately, hence the default ``B-E``.
    fontsize : float, optional
        Font size for all text.
    label_position : tuple, optional
        Axes-fraction position of the panel letters.
    title : str, optional
        Figure-level title (used for the bin-width series).
    ev_title : str, optional
        Title of the first panel.
    save_path : path-like, optional
        If given, the figure is saved here (PDF/PNG by extension).

    Returns
    -------
    (fig, axd)
        The figure and the dict of axes keyed by panel label.
    """
    if dihedral not in ("2D", "phi", "psi"):
        raise ValueError("dihedral must be '2D', 'phi' or 'psi'")

    labels = list(panel_labels)
    mosaic = f"{labels[0] * 2}{labels[1] * 2};{labels[2] * 2}{labels[3] * 2}"
    fig, axd = plt.subplot_mosaic(
        mosaic,
        sharex=False,
        figsize=figsize,
        dpi=dpi,
        facecolor="w",
        gridspec_kw={"height_ratios": [2, 2]},
    )

    hellingers = qual.hellingers_distances()
    if dihedral == "2D":
        metric = np.asarray(hellingers)
        all_to_all = np.asarray(qual.get_all_to_all_2d_trj_comparison())
    else:
        idx_sel = 0 if dihedral == "phi" else 1
        metric = np.asarray(hellingers[idx_sel])
        all_to_all = np.asarray(qual.get_all_to_all_trj_comparisons()[idx_sel])

    trj_helicity, _ = qual.fractional_helicity()
    trj_helicity = np.asarray(trj_helicity)

    # dihedral columns map to residues via residue_indices (0-based, caps
    # counted); express them 1-based without caps so capped and uncapped
    # chains are labelled the same way. Helicity (DSSP, caps excluded) is
    # always one column per non-cap residue.
    protein = qual.trajs[0].proteinTrajectoryList[qual.proteinID]
    offset = 1 if protein.ncap else 0
    idx = np.asarray(qual.residue_indices) - offset + 1
    helix_idx = np.arange(1, trj_helicity.shape[-1] + 1)
    xmax = max(idx[-1], helix_idx[-1])
    xticks = np.arange(increment, xmax + 1, increment)
    yticks = [0, 0.2, 0.4, 0.6, 0.8, 1]

    def _style(ax, ylabel, ax_title=None, xlabel=None):
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, fontsize=fontsize)
        ax.set_ylim([0, 1])
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=fontsize)
        ax.set_xlim([0, xmax + 1])
        if ax_title:
            ax.set_title(ax_title, fontsize=fontsize)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=fontsize)

    ax = axd[labels[0]]
    _style(ax, "Hellinger's Distance", ev_title)
    ax.plot(idx, metric.T, ".r", ms=4, alpha=0.3, mew=0)
    ax.plot(idx, metric.mean(axis=0), "sk-", ms=2, mew=0, linewidth=0.5)

    ax = axd[labels[1]]
    _style(ax, "Hellinger's Distance", "All-to-All Trajectory Comparison")
    ax.plot(idx, all_to_all.T, ".r", ms=4, alpha=0.3, mew=0)
    ax.plot(idx, all_to_all.mean(axis=0), "sk-", ms=2, mew=0, linewidth=0.5)

    ax = axd[labels[2]]
    _style(ax, "Hellinger's Distance\nmax - min", ev_title, xlabel="Residue")
    ax.bar(idx, np.ptp(metric, axis=0), width=0.8, color="k")

    ax = axd[labels[3]]
    _style(
        ax, "Fractional Helicity", "Secondary Structure Propensity", xlabel="Residue"
    )
    ax.plot(helix_idx, trj_helicity.T, ".r", ms=4, alpha=0.3, mew=0)
    ax.plot(helix_idx, trj_helicity.mean(axis=0), "sk-", ms=2, mew=0, linewidth=0.5)

    for label in labels:
        _panel_label(fig, axd[label], label, fontsize, label_position)
        if highlight is not None:
            axd[label].axvspan(
                highlight[0], highlight[1], alpha=0.5, facecolor=color, lw=0
            )
        if divider is not None:
            axd[label].axvline(divider, linewidth=0.5, linestyle="--", color=color)

    if title:
        fig.suptitle(title, fontsize=fontsize)

    fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi)
    return fig, axd


def print_sequence(qual: SamplingQuality, width: int = 10) -> str:
    """Return (and print) the simulated sequence with residue numbering.

    The manuscript draws the sequence as panel A of each figure in a vector
    graphics editor; here we simply print it with a ruler so the residue
    numbers in the plots can be mapped back to amino acids. Caps are removed.

    Parameters
    ----------
    qual : SamplingQuality
    width : int, optional
        Number of residues per ruler tick.

    Returns
    -------
    str
        The one-letter sequence without caps.
    """
    seq = (
        qual.trajs[0]
        .proteinTrajectoryList[qual.proteinID]
        .get_amino_acid_sequence(oneletter=True)
        .replace("<", "")
        .replace(">", "")
    )
    ruler = "".join(
        str((i + 1) // width % 10) if (i + 1) % width == 0 else "."
        for i in range(len(seq))
    )
    print(seq)
    print(ruler)
    return seq


def mean_bin_width_curves(
    build: "callable",
    bin_widths: Iterable[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the replicate-mean PENGUIN curves for a series of bin widths.

    Parameters
    ----------
    build : callable
        ``build(bwidth_degrees) -> SamplingQuality``.
    bin_widths : iterable of float
        Bin widths in degrees.

    Returns
    -------
    (bin_widths, ev_means, all_to_all_means)
        Arrays of shape ``(n_widths,)``, ``(n_widths, n_residues)`` and
        ``(n_widths, n_residues)``.
    """
    ev_means, ata_means = [], []
    widths = list(bin_widths)
    for bw in widths:
        qual = build(bw)
        ev_means.append(np.mean(np.asarray(qual.hellingers_distances()), axis=0))
        ata_means.append(
            np.mean(np.asarray(qual.get_all_to_all_2d_trj_comparison()), axis=0)
        )
    return np.array(widths), np.array(ev_means), np.array(ata_means)


def plot_bin_width_comparison(
    bin_widths: np.ndarray,
    curves: np.ndarray,
    title: str,
    save_path: Optional[os.PathLike] = None,
    figsize: Tuple[float, float] = (7, 5),
    fontsize: float = 10,
):
    """Overlay replicate-mean Hellinger curves for several bin widths.

    Parameters
    ----------
    bin_widths : array_like
        Bin widths in degrees, one per row of ``curves``.
    curves : np.ndarray
        Shape ``(n_widths, n_residues)``.
    title : str
        Axes title.
    save_path : path-like, optional
        Where to save the figure.
    figsize, fontsize : optional

    Returns
    -------
    (fig, ax)
    """
    residues = np.arange(1, curves.shape[1] + 1)
    fig, ax = plt.subplots(figsize=figsize, facecolor="w")
    for bw, curve in zip(bin_widths, curves):
        ax.plot(residues, curve, label=f"Bin width: {bw:.1f}°")
    ax.set_ylim([0, 1])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel("Residue", fontsize=fontsize)
    ax.set_ylabel("Hellinger's Distance", fontsize=fontsize)
    ax.legend(frameon=False, fontsize=fontsize)
    fig.tight_layout()
    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300)
    return fig, ax


# ---------------------------------------------------------------------------
# Synthetic energy landscapes (Figures 1 and 2B)
# ---------------------------------------------------------------------------


def _gaussian_multD(xs, x0s, height=1.0, widths=1.0):
    """Unnormalised multi-dimensional Gaussian evaluated on a grid."""
    xs = np.asarray(xs)
    x0s = np.broadcast_to(np.asarray(x0s, dtype=float), (len(xs),))
    widths = np.broadcast_to(np.asarray(widths, dtype=float), (len(xs),))
    exponent = np.sum(
        [((xs[k] - x0s[k]) ** 2) / (2 * widths[k] ** 2) for k in range(len(xs))],
        axis=0,
    )
    return height * np.exp(-exponent)


def gaussian_noise(
    x1s: np.ndarray,
    x2s: np.ndarray,
    gaussians_per_axis=None,
    height_range=(-0.1, 0.1),
    width_range=(0.9, 1.1),
    rigidity: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Sum of randomly placed 2D Gaussians used to draw cartoon energy landscapes.

    Adapted from ``slandscapes`` (M. Zimmerman,
    https://github.com/mizimmer90/slandscapes). A grid of
    ``gaussians_per_axis`` x ``gaussians_per_axis`` Gaussian bumps is laid over
    the ``(x1s, x2s)`` mesh; each bump gets a random height and width from the
    given ranges and is jittered off its grid position by an amount controlled
    by ``rigidity`` (0 = loose, 1 = rigid).

    Parameters
    ----------
    x1s, x2s : np.ndarray
        Output of ``np.meshgrid``.
    gaussians_per_axis : int or (int, int), optional
        Number of bumps along each axis. Default: half the axis length.
    height_range, width_range : tuple of float, optional
        Uniform sampling ranges for bump height and width.
    rigidity : float, optional
        How tightly bumps stick to their grid positions.
    rng : np.random.Generator, optional
        Random generator; ``np.random.default_rng()`` if omitted. The
        manuscript figures used the legacy ``np.random.seed`` API; pass an
        explicit generator for reproducibility.

    Returns
    -------
    np.ndarray
        Landscape values on the mesh.
    """
    rng = np.random.default_rng() if rng is None else rng
    if isinstance(gaussians_per_axis, int):
        gaussians_per_axis = [gaussians_per_axis, gaussians_per_axis]
    elif gaussians_per_axis is None:
        gaussians_per_axis = [int(x1s.max() / 2.0), int(x2s.max() / 2.0)]
    tot = int(np.prod(gaussians_per_axis))

    axis1 = x1s[0]
    axis2 = x2s[:, 0]
    box1 = (axis1[-1] - axis1[0]) / gaussians_per_axis[0]
    box2 = (axis2[-1] - axis2[0]) / gaussians_per_axis[1]
    centers1 = np.arange(gaussians_per_axis[0]) * box1 + box1 / 2.0
    centers2 = np.arange(gaussians_per_axis[1]) * box2 + box2 / 2.0
    centers = np.array([(a, b) for a in centers1 for b in centers2])

    heights = rng.uniform(height_range[0], height_range[1], size=tot)
    widths = rng.uniform(width_range[0], width_range[1], size=tot)
    rigidity_div = 2 + (rigidity * 4) ** 2

    noise = np.zeros((len(axis2), len(axis1)))
    xs = np.array([x1s, x2s])
    for k in range(tot):
        jitter = rng.random(2) - 0.5
        center = [
            centers[k, 0] + (box1 / rigidity_div) * jitter[0],
            centers[k, 1] + (box2 / rigidity_div) * jitter[1],
        ]
        noise += _gaussian_multD(xs, center, height=heights[k], widths=widths[k])
    return noise
