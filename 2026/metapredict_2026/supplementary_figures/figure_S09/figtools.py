"""Shared figure style, colours and paths for the metapredict V3 manuscript figures.

Every figure directory carries an identical copy of this file so that each one runs
on its own, while every panel in the manuscript draws metapredict V1/V2/V3 in the
same colours and with the same text settings. The master copy lives in
``shared/figtools.py``; edit that and re-copy it with ``shared/sync_figtools.sh``
rather than editing one directory's copy.

Figure directories hold code only. Every input a notebook reads lives in ``data/`` at
the top of the repository and is found with :func:`shared_data_path` (benchmark
results with :func:`benchmark_results_path`).
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------- #
# Colours
# ----------------------------------------------------------------------------- #
# metapredict version colours, used for V1/V2/V3 in every panel. These are the
# hues from the original CAID ranking notebook (chosen so V1/V2/V3 are easy to
# tell apart).
V1_COLOR = "#FF9D33"  # orange
V2_COLOR = "#2978d6"  # blue
V3_COLOR = "#CC2D46"  # crimson
VERSION_COLORS = {"V1": V1_COLOR, "V2": V2_COLOR, "V3": V3_COLOR}

# Default categorical order (xmgrace) for up to five series that are NOT
# metapredict versions (e.g. hardware types, sequence sets).
XMGRACE_COLORS = ["black", "red", "blue", "green", "purple"]

# Bars for "the rest of the field" in ranking panels.
FIELD_COLOR = "black"

# Fixed colours for the top CAID3 Disorder-PDB predictors whenever they are drawn
# next to metapredict (Figure 2, Figure S4), so each method keeps one colour across
# panels. The xmgrace red/blue are deliberately avoided here because they read as
# metapredict V3/V2.
CAID3_PREDICTOR_COLORS = {
    "PUNCH2": "#1baf7a",          # green
    "PUNCH2-Light": "#6a3d9a",    # violet
    "AlphaFold-rsa": "#0e8f9c",   # teal
    "SPOT-Disorder2": "#111111",  # near-black
    "AlphaFold3-rsa": "#7a4d1e",  # brown
}

# CAID3 Disorder-PDB reference labels, drawn as a thin bar under per-residue traces
# (Figure 2C-E, S4, S10). Keys are the reference label values: 1 disordered,
# 0 ordered, -1 unscored (drawn hatched).
DISORDER_PDB_COLORS = {
    1: "#e8853a",   # orange
    0: "#1baf7a",   # green
    -1: "#d7d7d3",  # light grey
}
HATCH_EDGE_COLOR = "#b4b4b0"  # hatch lines on "no information" blocks

# ----------------------------------------------------------------------------- #
# Paths
# ----------------------------------------------------------------------------- #
# The repository root is two levels above every figure directory
# (<repo>/<maintext|supplementary>_figures/<figure>/figtools.py).
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Every input the figures read (CAID predictions and reference sets, PDB coverage,
# sequences, AlphaFold models and renders, viral-proteome tables, benchmark
# results); see data/README.md. Override with an environment variable if the data
# live elsewhere.
DATA_ENV_VAR = "METAPREDICT_V3_DATA_DIR"
SHARED_DATA_DIR = Path(os.environ.get(DATA_ENV_VAR, REPO_ROOT / "data"))

# TSVs written by the benchmark scripts in benchmarks/ (run those on the benchmark
# machines), one subdirectory per benchmark.
BENCHMARK_RESULTS_DIR = SHARED_DATA_DIR / "benchmarks"

# Panel PDFs are written here, relative to the notebook's working directory.
FIGURE_DIR = Path("figures")

# Bold panel letters (A, B, ...) in the combined multi-panel figures, in points.
PANEL_LETTER_FONTSIZE = 12


def shared_data_path(*parts: str) -> Path:
    """Return a path inside the shared data directory, failing loudly if absent.

    Parameters
    ----------
    *parts : str
        Path components below ``data/`` (or ``$METAPREDICT_V3_DATA_DIR``), e.g.
        ``("caid_references", "caid3_disorder_pdb.fasta")``.

    Returns
    -------
    pathlib.Path
        The existing path.

    Raises
    ------
    FileNotFoundError
        If the path does not exist, naming what was expected and where to get it.
    """
    path = SHARED_DATA_DIR.joinpath(*parts)
    if not path.exists():
        raise FileNotFoundError(
            f"MISSING: expected {path}. See data/README.md for how to "
            f"download or regenerate it (or set ${DATA_ENV_VAR})."
        )
    return path


def metapredict_caid_path(caid_round: str, version: str) -> Path:
    """Path to the metapredict 3.1.0 predictions for one CAID round and network.

    These are built by ``data/build_metapredict_caid_predictions.py`` and replace
    the metapredict files distributed by CAID (made with earlier versions), so every
    metapredict number in the manuscript comes from 3.1.0. Binary calls use the
    3.1.0 default (residues inside predicted IDRs are disordered).

    Parameters
    ----------
    caid_round : str
        ``"caid1"``, ``"caid2"`` or ``"caid3"``.
    version : str
        ``"V1"``, ``"V2"`` or ``"V3"``.

    Returns
    -------
    pathlib.Path
        Path to ``Metapredict-v<N>.caid``.
    """
    return shared_data_path("metapredict_caid", "domains", caid_round,
                            f"Metapredict-v{version[-1]}.caid")


def benchmark_results_path(benchmark: str) -> Path:
    """Return the results directory of one benchmark, failing loudly if absent.

    Parameters
    ----------
    benchmark : str
        Benchmark name, which is also the subdirectory of ``benchmarks/`` holding its
        script: ``"proteome"``, ``"batch_size"``, ``"batch_vs_sequential"`` or
        ``"stream"``.

    Returns
    -------
    pathlib.Path
        ``data/benchmarks/<benchmark>``, the directory that benchmark's script writes
        its TSVs to.

    Raises
    ------
    FileNotFoundError
        If the directory does not exist.
    """
    path = BENCHMARK_RESULTS_DIR / benchmark
    if not path.exists():
        raise FileNotFoundError(
            f"MISSING: expected benchmark results in {path}. Run the script in "
            f"benchmarks/{benchmark}/ first (see benchmarks/README.md)."
        )
    return path


# ----------------------------------------------------------------------------- #
# metapredict version guard
# ----------------------------------------------------------------------------- #
# Every figure is made with metapredict 3.1.x.
REQUIRED_METAPREDICT_SERIES = "3.1"


def require_metapredict_3_1() -> str:
    """Check that metapredict 3.1.x is the version being imported, and report it.

    Prints the version, the path it is imported from and, when that path is a git
    checkout, the commit, so each executed notebook records exactly which code made
    its predictions.

    Returns
    -------
    str
        The metapredict version string.

    Raises
    ------
    RuntimeError
        If the imported metapredict is not a 3.1.x release.
    """
    import subprocess

    import metapredict

    version = metapredict.__version__
    location = Path(metapredict.__file__).resolve().parent
    if not version.startswith(REQUIRED_METAPREDICT_SERIES):
        raise RuntimeError(
            f"metapredict {version} found at {location}; these figures require "
            f"metapredict {REQUIRED_METAPREDICT_SERIES}.x"
        )
    commit = subprocess.run(
        ["git", "-C", str(location), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True,
    ).stdout.strip()
    print(f"metapredict {version} from {location}" + (f" (git {commit})" if commit else ""))
    return version


# ----------------------------------------------------------------------------- #
# Style
# ----------------------------------------------------------------------------- #
def apply_style() -> None:
    """Apply the manuscript matplotlib style.

    Arial throughout, text kept as editable text in PDF/PS/SVG exports (so it can
    be edited in Illustrator), thin solid-black axes, and no grid lines.
    """
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "axes.linewidth": 0.5,
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "axes.grid": False,
            "xtick.color": "black",
            "ytick.color": "black",
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "legend.frameon": False,
            "savefig.bbox": "tight",
        }
    )


def save_panel(fig: plt.Figure, name: str) -> Path:
    """Save a panel as ``figures/<name>.pdf`` and return the path.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    name : str
        File stem, e.g. ``"figure_1A"``.

    Returns
    -------
    pathlib.Path
        Path of the written PDF.
    """
    FIGURE_DIR.mkdir(exist_ok=True)
    out = FIGURE_DIR / f"{name}.pdf"
    fig.savefig(out)
    print(f"Wrote {out}")
    return out


def label_panels(fig: plt.Figure, axes: list, letters: str) -> None:
    """Write a bold panel letter at the top-left corner of each panel in a combined figure.

    Each letter is placed at the top-left of everything drawn for that axes (tick
    labels, axis labels and any text above it), so it never covers the panel. Call
    this after the layout is final (after ``tight_layout``, or on a constrained-layout
    figure), because it reads where the axes have been drawn.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The combined figure; letters are added as figure text.
    axes : list of matplotlib.axes.Axes
        One axes per panel, in letter order (for a panel made of several axes, pass
        its top axes).
    letters : str
        The panel letters, one per axes, e.g. ``"ABCD"``.
    """
    assert len(axes) == len(letters), f"{len(axes)} axes but {len(letters)} letters"
    fig.draw_without_rendering()
    display_to_figure = fig.transFigure.inverted()
    for ax, letter in zip(axes, letters):
        extent = ax.get_tightbbox().transformed(display_to_figure)
        fig.text(extent.x0, extent.y1, letter, fontsize=PANEL_LETTER_FONTSIZE, fontweight="bold",
                 ha="left", va="top")
