"""Loading and plotting helpers for the metapredict performance figures.

Identical copies of this file live in maintext_figures/figure_1 (panel 1E) and
supplementary_figures/figure_S05, figure_S06, figure_S07 and figure_S08. Each
reads the TSVs that the scripts in ``benchmarks/`` write to ``data/benchmarks/`` and plots
whatever (version, machine, device) combinations are present, so results from new
machines drop in without code changes.

Rows from a metapredict other than 3.1.x are dropped with a warning, so versions
are never silently mixed in one panel.
"""

from __future__ import annotations

import csv
import re
import warnings
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter

from figtools import (
    REQUIRED_METAPREDICT_SERIES,
    VERSION_COLORS,
    XMGRACE_COLORS,
    benchmark_results_path,
)

VERSION_ORDER = ["V1", "V2", "V3"]

# Hardware groups are ordered CPUs first, then GPUs.
DEVICE_ORDER = {"cpu": 0, "mps": 1, "cuda": 2}
DEVICE_NAMES = {"cpu": "CPU", "mps": "GPU (MPS)", "cuda": "GPU (CUDA)"}

# Non-version series take the xmgrace colour order.
DEVICE_COLORS = {"cpu": XMGRACE_COLORS[0], "mps": XMGRACE_COLORS[1], "cuda": XMGRACE_COLORS[2]}

# Log-scale time axes are labelled at 1, 2, 5 × each decade with plain numbers.
LOG_TICK_SUBS = (1.0, 2.0, 5.0)


# ----------------------------------------------------------------------------- #
# Loading
# ----------------------------------------------------------------------------- #
def load_results(benchmark: str) -> list[dict]:
    """Read every results TSV of one benchmark, keeping only metapredict 3.1.x rows.

    Parameters
    ----------
    benchmark : str
        Benchmark name (``proteome``, ``batch_size``, ``batch_vs_sequential`` or
        ``stream``), i.e. the subdirectory of ``data/benchmarks`` holding its TSVs.

    Returns
    -------
    list of dict
        One dict per row (values are strings, as read); ``source_file`` is added.

    Raises
    ------
    FileNotFoundError
        If there are no TSVs, or no rows survive the version filter.
    """
    files = sorted(benchmark_results_path(benchmark).glob("*.tsv"))
    if not files:
        raise FileNotFoundError(f"MISSING: no TSV files in {benchmark_results_path(benchmark)}")
    rows = []
    for path in files:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row["source_file"] = path.name
                rows.append(row)

    kept = [r for r in rows if r["metapredict_version"].startswith(REQUIRED_METAPREDICT_SERIES)]
    dropped = sorted({(r["source_file"], r["metapredict_version"]) for r in rows} -
                     {(r["source_file"], r["metapredict_version"]) for r in kept})
    if dropped:
        warnings.warn(f"dropped {len(rows) - len(kept)} rows not from metapredict "
                      f"{REQUIRED_METAPREDICT_SERIES}.x: {dropped}")
    if not kept:
        raise FileNotFoundError(f"no metapredict {REQUIRED_METAPREDICT_SERIES}.x rows in {benchmark}")

    # Smoke-test files (tiny inputs) are only used when there is nothing else.
    smoke = sorted({r["source_file"] for r in kept if r["source_file"].startswith("SMOKETEST_")})
    real = [r for r in kept if not r["source_file"].startswith("SMOKETEST_")]
    if smoke and real:
        print(f"  ignoring {len(smoke)} SMOKETEST_ files because real results exist")
        kept = real
    elif smoke:
        warnings.warn(f"plotting SMOKE-TEST data ({len(smoke)} files); replace with real runs")
    versions = sorted({(r["metapredict_version"], r["metapredict_git_commit"]) for r in kept})
    print(f"{benchmark}: {len(kept)} rows from {len({r['source_file'] for r in kept})} files; "
          f"metapredict (version, commit): {versions}")
    return kept


def values(rows: list[dict], key: str) -> np.ndarray:
    """Return column ``key`` of ``rows`` as a float array (empty strings -> NaN)."""
    return np.array([float(r[key]) if r[key] != "" else np.nan for r in rows])


def _clean_name(name: str) -> str:
    """Shorten a CPU/GPU brand string, e.g. 'Apple M3 Max' -> 'M3 Max'."""
    for junk in ["(R)", "(TM)", "Apple ", "NVIDIA ", "GeForce ", "CPU", "Processor"]:
        name = name.replace(junk, "")
    name = re.sub(r"@.*", "", name)
    name = re.sub(r"\d+-Core", "", name)
    return " ".join(name.split())


def machine_name(row: dict) -> str:
    """Name of the machine that produced a row: its label, else its CPU brand."""
    return row["machine_label"] or _clean_name(row["cpu_brand"]) or row["hostname"]


def hardware_label(row: dict) -> str:
    """Two-line hardware label for a row, e.g. 'M3 Max\\nGPU (MPS)'.

    CUDA rows are named by their GPU; CPU and MPS rows by their machine.
    """
    device = row["device_used"]
    name = _clean_name(row["gpu_name"]) if device == "cuda" else machine_name(row)
    return f"{name}\n{DEVICE_NAMES[device]}"


def hardware_groups(rows: list[dict]) -> list[str]:
    """Distinct hardware labels in rows, CPUs first then MPS then CUDA."""
    seen = {hardware_label(r): DEVICE_ORDER[r["device_used"]] for r in rows}
    return sorted(seen, key=lambda label: (seen[label], label))


def median_by(rows: list[dict], key_fn, value_key: str) -> dict:
    """Median of ``value_key`` for each group defined by ``key_fn(row)``.

    Repeated runs of the same configuration are reported and summarised by their
    median rather than silently overwritten.

    Parameters
    ----------
    rows : list of dict
        Benchmark rows.
    key_fn : callable
        Maps a row to a hashable group key.
    value_key : str
        Column to summarise.

    Returns
    -------
    dict
        ``{group key: median value}``; rows with an empty value are skipped.
    """
    groups = defaultdict(list)
    for row in rows:
        if row[value_key] != "":
            groups[key_fn(row)].append(float(row[value_key]))
    for key, vals in groups.items():
        if len(vals) > 1:
            print(f"  {key}: {len(vals)} repeat runs, using the median")
    return {key: float(np.median(vals)) for key, vals in groups.items()}


# ----------------------------------------------------------------------------- #
# Proteome benchmark (Figure 1E, Figure S5)
# ----------------------------------------------------------------------------- #
def proteome_rows() -> list[dict]:
    """Load the proteome benchmark, keeping score-only runs over the full proteome.

    Runs with ``--limit`` (sequence subsets) or ``--domain-decomposition`` are
    dropped when full score-only runs exist; if none do, the limited runs are kept
    with a warning (smoke-test data).
    """
    rows = load_results("proteome")
    rows = [r for r in rows if r["domain_decomposition"] == "False"]
    full = [r for r in rows if r["sequence_limit"] == ""]
    if full:
        if len(full) < len(rows):
            print(f"  dropped {len(rows) - len(full)} runs on sequence subsets (--limit)")
        return full
    warnings.warn("no full-proteome runs found; plotting runs on sequence subsets")
    return rows


def version_legend(ax: plt.Axes, **kwargs) -> None:
    """Add a V1/V2/V3 colour legend to ``ax``."""
    handles = [Patch(facecolor=VERSION_COLORS[v], edgecolor="black", linewidth=0.5, label=v)
               for v in VERSION_ORDER]
    ax.legend(handles=handles, fontsize=7, handlelength=1.2, columnspacing=1.0, **kwargs)


def grouped_version_bars(ax: plt.Axes, groups: list[str], value_by_key: dict,
                         label_fmt=lambda v: f"{v:.0f}", label_rotation: int = 0,
                         group_gap: float = 1.0, bottom: float = 0.0) -> np.ndarray:
    """Draw V1/V2/V3 bars side by side for each hardware group.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    groups : list of str
        Group labels, left to right.
    value_by_key : dict
        ``{(version, group label): value}``; missing combinations are left empty.
    label_fmt : callable, optional
        Formats the value written above each bar.
    label_rotation : int, optional
        Rotation of the bar labels in degrees.
    group_gap : float, optional
        Extra space between groups, in bar widths.
    bottom : float, optional
        Value the bars start from. On a log axis pass the lower y-limit: a bar from 0
        starts at minus infinity, which matplotlib writes into the PDF as a rectangle
        reaching tens of thousands of points below the page (clipped on screen, but
        Illustrator refuses to paste it).

    Returns
    -------
    numpy.ndarray
        x positions of the group centres (for tick labels).
    """
    centres = []
    for g_idx, group in enumerate(groups):
        base = g_idx * (len(VERSION_ORDER) + group_gap)
        xs = base + np.arange(len(VERSION_ORDER))
        centres.append(xs.mean())
        for x, version in zip(xs, VERSION_ORDER):
            value = value_by_key.get((version, group))
            if value is None:
                continue
            assert value > bottom, f"{version} {group}: value {value} is not above the bar bottom {bottom}"
            ax.bar(x, value - bottom, 0.8, bottom=bottom, color=VERSION_COLORS[version],
                   edgecolor="black", linewidth=0.5, zorder=3)
            ax.annotate(label_fmt(value), (x, value), textcoords="offset points",
                        xytext=(0, 2), ha="center", va="bottom", fontsize=6.5,
                        rotation=label_rotation)
    ax.set_xticks(centres, groups, fontsize=7)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=7)
    return np.array(centres)


def gpu_speedups(rows: list[dict]) -> dict:
    """GPU-over-CPU throughput ratio per (version, machine) with both devices.

    Parameters
    ----------
    rows : list of dict
        Proteome benchmark rows.

    Returns
    -------
    dict
        ``{(version, 'machine\\nMPS vs CPU'): speed-up}``, pairing each GPU run with
        the CPU run of the same version on the same host.
    """
    throughput = median_by(rows, lambda r: (r["network_version"], r["hostname"], r["device_used"]),
                           "residues_per_second")
    names = {r["hostname"]: (machine_name(r), _clean_name(r["gpu_name"])) for r in rows if r["gpu_name"]}
    speedups = {}
    for (version, host, device), rate in throughput.items():
        cpu_rate = throughput.get((version, host, "cpu"))
        if device == "cpu" or cpu_rate is None:
            continue
        machine, gpu = names[host]
        name = gpu if device == "cuda" else machine
        speedups[(version, f"{name}\n{device.upper()} vs CPU")] = rate / cpu_rate
    return speedups


# ----------------------------------------------------------------------------- #
# Batch-size sweep (Figure S6)
# ----------------------------------------------------------------------------- #
def plot_batch_size_panel(ax: plt.Axes, rows: list[dict], version: str,
                          show_legend: bool = True) -> None:
    """Prediction time vs batch size for one network, one line per machine/device.

    The batch size metapredict 3.1 uses by default on each device is marked with
    an open square. Batch sizes that ran out of GPU memory have no row in the
    results, so their line simply stops.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    rows : list of dict
        Batch-size sweep rows (all versions).
    version : str
        ``"V1"``, ``"V2"`` or ``"V3"``.
    show_legend : bool, optional
        Draw the legend. The series are the same in every panel, so the assembled
        figure only needs it once.
    """
    rows = [r for r in rows if r["network_version"] == version]
    assert rows, f"no batch-size rows for {version}"
    # only compare sweeps over the same input (the largest one, normally the full proteome)
    n_max = max(int(r["n_sequences"]) for r in rows)
    if any(int(r["n_sequences"]) != n_max for r in rows):
        print(f"  dropped {sum(int(r['n_sequences']) != n_max for r in rows)} rows run on "
              f"fewer than {n_max} sequences")
        rows = [r for r in rows if int(r["n_sequences"]) == n_max]
    times = median_by(rows, lambda r: (hardware_label(r), r["device_used"], int(r["batch_size"])),
                      "total_seconds")
    defaults = {(hardware_label(r), r["device_used"]): int(r["default_batch_size"]) for r in rows}

    linestyles = ["-", "--", ":", "-."]
    machine_index = {m: i for i, m in enumerate(sorted({machine_name(r) for r in rows}))}
    for (label, device), default_bs in sorted(defaults.items(), key=lambda kv: DEVICE_ORDER[kv[0][1]]):
        sizes = sorted(bs for (lab, dev, bs) in times if (lab, dev) == (label, device))
        secs = [times[(label, device, bs)] for bs in sizes]
        machine = next(machine_name(r) for r in rows if hardware_label(r) == label)
        ax.plot(sizes, secs, linestyles[machine_index[machine] % len(linestyles)], marker="o",
                markersize=3, color=DEVICE_COLORS[device], linewidth=1, label=label.replace("\n", " "))
        if default_bs in sizes:
            ax.plot(default_bs, times[(label, device, default_bs)], "s", markersize=8,
                    markerfacecolor="none", markeredgecolor=DEVICE_COLORS[device], markeredgewidth=1)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    # plain-number ticks (5, 10, 20, ...) instead of mixed 10^n / 2×10^n labels
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=LOG_TICK_SUBS))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([32, 64, 128, 256, 512, 1024], ["32", "64", "128", "256", "512", "1024"])
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Batch size (sequences)", fontsize=9)
    ax.set_ylabel(f"{version} prediction time (s)", fontsize=9)
    if show_legend:
        ax.plot([], [], "s", markersize=6, markerfacecolor="none", markeredgecolor="black",
                label="3.1.0 default")
        ax.legend(fontsize=6.5, loc="upper right")


# ----------------------------------------------------------------------------- #
# Batch vs sequential (Figure S7)
# ----------------------------------------------------------------------------- #
# Sequential predictions (always on the CPU) are red; batched predictions are
# black on a CPU and blue on a GPU.
SEQUENTIAL_COLOR = XMGRACE_COLORS[1]
BATCH_COLORS = {"cpu": XMGRACE_COLORS[0], "mps": XMGRACE_COLORS[2], "cuda": XMGRACE_COLORS[2]}


def plot_batch_vs_sequential_panel(ax: plt.Axes, rows: list[dict], version: str) -> None:
    """Time vs number of sequences for batched and one-at-a-time prediction.

    Points are measurements (median over replicates); dashed lines are linear fits,
    as in the original figure. One marker shape per machine. A machine with a GPU
    can contribute batched runs on both its CPU and its GPU; these are separate
    series, never pooled.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    rows : list of dict
        Batch-vs-sequential rows (all versions).
    version : str
        ``"V1"``, ``"V2"`` or ``"V3"``.
    """
    rows = [r for r in rows if r["network_version"] == version]
    assert rows, f"no batch-vs-sequential rows for {version}"

    def series_key(row: dict) -> tuple[str, str, str]:
        return machine_name(row), row["mode"], hardware_label(row)

    times = median_by(rows, lambda r: (*series_key(r), int(r["n_sequences"])), "total_seconds")
    machines = sorted({machine_name(r) for r in rows})
    devices = {hardware_label(r): r["device_used"] for r in rows}
    # per machine: sequential first, then batched CPU, then batched GPU
    series = sorted({key[:3] for key in times},
                    key=lambda s: (machines.index(s[0]), s[1] != "sequential", DEVICE_ORDER[devices[s[2]]]))

    markers = ["o", "^", "s", "D"]
    for machine, mode, hardware in series:
        counts = sorted(n for (mach, md, hw, n) in times if (mach, md, hw) == (machine, mode, hardware))
        secs = [times[(machine, mode, hardware, n)] for n in counts]
        color = SEQUENTIAL_COLOR if mode == "sequential" else BATCH_COLORS[devices[hardware]]
        hardware_name = hardware.replace("\n", " ")
        ax.plot(counts, secs, markers[machines.index(machine) % len(markers)], color=color, ms=4,
                label=f"{mode}, {hardware_name}")
        slope, intercept = np.polyfit(counts, secs, 1)
        fit_x = np.array([0, max(counts)])
        ax.plot(fit_x, slope * fit_x + intercept, "--", color=color, lw=0.5)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Number of sequences", fontsize=8)
    ax.set_ylabel(f"{version} prediction time (s)", fontsize=8)
    ax.legend(fontsize=6, loc="upper left")


# ----------------------------------------------------------------------------- #
# Streaming memory (Figure S8)
# ----------------------------------------------------------------------------- #
STREAM_COLORS = {"full": XMGRACE_COLORS[0], "stream": XMGRACE_COLORS[1]}
STREAM_LABELS = {"full": "predict_disorder", "stream": "predict_disorder_stream"}

# Human reference proteome size (UP000005640, 20,652 sequences), for a reference line.
HUMAN_PROTEOME_RESIDUES = 11_453_425

# Machine-RAM ceilings the extrapolated full-mode curve is projected onto (GiB).
RAM_LINES_GB = [(16.0, "16 GB laptop"), (64.0, "64 GB workstation")]


# Measured stream series are solid lines with one marker shape per machine and
# device (dashed lines are reserved for the projection).
STREAM_MARKERS = ["o", "^", "s", "D"]


def stream_hardware(sweep: list[dict]) -> str:
    """One-line hardware name of a stream sweep, e.g. 'RTX A4000 GPU (CUDA)'."""
    return hardware_label(sweep[0]).replace("\n", " ")


def stream_sweeps(rows: list[dict]) -> list[list[dict]]:
    """Pick one stream sweep per machine and device, the most complete sweep first.

    For each machine and device, prefers the full sweep, then the quick sweep, then
    custom (smoke-test) sizes, and within that the most recent file. Sweeps are
    ordered by completeness, then CPUs before GPUs, then by name; the RAM
    projection is fitted to the first one alone. Prints which files were chosen.

    Parameters
    ----------
    rows : list of dict
        Stream benchmark rows (all machines and devices).

    Returns
    -------
    list of list of dict
        One list of rows per machine and device.
    """
    preference = {"full": 0, "quick": 1, "custom": 2}
    files_by_hardware = defaultdict(set)
    sweep_rank = {}
    device_rank = {}
    latest = defaultdict(str)
    for r in rows:
        files_by_hardware[hardware_label(r)].add(r["source_file"])
        sweep_rank[r["source_file"]] = preference[r["sweep"]]
        device_rank[r["source_file"]] = (DEVICE_ORDER[r["device_used"]], hardware_label(r))
        latest[r["source_file"]] = max(latest[r["source_file"]], r["timestamp_utc"])

    chosen = []
    for files in files_by_hardware.values():
        best = min(sweep_rank[f] for f in files)
        chosen.append(max((f for f in files if sweep_rank[f] == best), key=latest.get))
    chosen.sort(key=lambda f: (sweep_rank[f], device_rank[f]))

    sweeps = [[r for r in rows if r["source_file"] == f] for f in chosen]
    for f, sweep in zip(chosen, sweeps):
        print(f"plotting stream sweep {f} ({stream_hardware(sweep)}, {len(sweep)} measurements)")
    return sweeps


def _stream_legend(ax: plt.Axes, sweeps: list[list[dict]], projected: bool, **legend_kwargs) -> None:
    """Legend that keys the function by colour and the hardware by marker shape.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to put the legend on.
    sweeps : list of list of dict
        The plotted sweeps, in plotting order (marker shapes follow this order).
    projected : bool
        Whether to add an entry for the dashed full-mode projection.
    **legend_kwargs
        Passed on to ``ax.legend`` (placement, font size, frame).
    """
    handles = [Line2D([], [], color=STREAM_COLORS[mode], lw=1.2, label=STREAM_LABELS[mode])
               for mode in ["full", "stream"]]
    if projected:
        handles.append(Line2D([], [], color=STREAM_COLORS["full"], lw=1, linestyle="--",
                              label=f"{STREAM_LABELS['full']}, projected"))
    for m_idx, sweep in enumerate(sweeps):
        handles.append(Line2D([], [], color="black", linestyle="none", markersize=4,
                              marker=STREAM_MARKERS[m_idx % len(STREAM_MARKERS)],
                              label=stream_hardware(sweep)))
    ax.legend(handles=handles, **legend_kwargs)


def _stream_series(rows: list[dict], mode: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(millions of residues, peak RSS in GB, throughput in 10^3 res/s) for one mode, sorted."""
    sub = sorted((r for r in rows if r["mode"] == mode), key=lambda r: int(r["n_residues"]))
    assert sub, f"no {mode} rows in the stream sweep"
    return (values(sub, "n_residues") / 1e6, values(sub, "peak_rss_gb"),
            values(sub, "residues_per_second") / 1e3)


def plot_stream_memory(ax: plt.Axes, sweeps: list[list[dict]]) -> None:
    """Peak RAM vs dataset size (log-log) with the full-mode projection to RAM limits.

    Full mode holds the whole input and output, so its peak grows linearly with the
    dataset; a least-squares line through the measured full-mode points of the
    first sweep is extrapolated (dashed) to where it crosses each RAM ceiling.
    Streaming is drawn as measured points only, for every sweep. On a GPU this is
    host RAM; GPU memory is not included. Logic from
    ``stream_performance/generate_stream_figure.py``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    sweeps : list of list of dict
        One stream sweep per machine and device, most complete first (see
        :func:`stream_sweeps`).
    """
    rows = sweeps[0]
    fx, fy, _ = _stream_series(rows, "full")
    sx, sy, _ = _stream_series(rows, "stream")
    slope, floor = np.polyfit(fx, fy, 1)  # GB per million residues, GB
    print(f"full-mode projection fitted to {stream_hardware(rows)}: {slope * 1e3:.2f} GB per billion "
          f"residues, floor {floor:.2f} GB")
    top_ram = max(gb for gb, _ in RAM_LINES_GB)

    ax.set_xscale("log")
    ax.set_yscale("log")
    if slope > 0:
        x_oom = (top_ram - floor) / slope  # dataset size where full mode exceeds top RAM
        x_hi = x_oom * 1.6
        # region only streaming can reach
        ax.axvspan(x_oom, x_hi, color=STREAM_COLORS["stream"], alpha=0.06, lw=0, zorder=0)
        for gb, label in RAM_LINES_GB:
            ax.axhline(gb, color="black", linewidth=0.5, linestyle=":", zorder=1)
            ax.annotate(label, (x_hi, gb), textcoords="offset points", xytext=(-2, 2),
                        ha="right", va="bottom", fontsize=6.5)
            ax.plot((gb - floor) / slope, gb, "v", color=STREAM_COLORS["full"], markersize=5, zorder=6)
        xs_proj = np.logspace(np.log10(fx.max()), np.log10(x_hi), 100)
        ax.plot(xs_proj, floor + slope * xs_proj, "--", color=STREAM_COLORS["full"], lw=1)
        ax.set_ylim(min(sy.min(), floor) * 0.7, top_ram * 2.0)
    else:
        warnings.warn("full-mode peak RAM does not grow with dataset size in this sweep; "
                      "no projection drawn (too few / too small datasets)")
        x_hi = max(fx.max(), sx.max()) * 1.6

    x_lo = min(fx.min(), sx.min())
    for m_idx, sweep in enumerate(sweeps):
        for mode in ["full", "stream"]:
            x, y, _ = _stream_series(sweep, mode)
            x_lo = min(x_lo, x.min())
            ax.plot(x, y, "-", color=STREAM_COLORS[mode], marker=STREAM_MARKERS[m_idx % len(STREAM_MARKERS)],
                    lw=1.2, markersize=4)

    hp = HUMAN_PROTEOME_RESIDUES / 1e6
    ax.axvline(hp, color="black", linewidth=0.5, zorder=1)
    ax.annotate("human\nproteome", (hp, ax.get_ylim()[1]), textcoords="offset points",
                xytext=(-3, -3), ha="right", va="top", fontsize=6.5)
    ax.set_xlim(min(x_lo, hp) * 0.6, x_hi)
    ax.set_xlabel("Dataset size (millions of residues)", fontsize=9)
    ax.set_ylabel("Peak RAM (GB)", fontsize=9)
    ax.tick_params(labelsize=7)
    # opaque, borderless legend so the human-proteome line does not run through it
    _stream_legend(ax, sweeps, projected=slope > 0, fontsize=6.3, loc="center left",
                   bbox_to_anchor=(0.02, 0.62), frameon=True, framealpha=1, edgecolor="none")


def plot_stream_throughput(ax: plt.Axes, sweeps: list[list[dict]]) -> None:
    """Throughput vs dataset size for predict_disorder and predict_disorder_stream.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    sweeps : list of list of dict
        One stream sweep per machine and device (see :func:`stream_sweeps`).
    """
    top = 0.0
    for m_idx, sweep in enumerate(sweeps):
        for mode in ["full", "stream"]:
            x, _, rate = _stream_series(sweep, mode)
            ax.plot(x, rate, "-", color=STREAM_COLORS[mode], marker=STREAM_MARKERS[m_idx % len(STREAM_MARKERS)],
                    lw=1.2, markersize=4)
            top = max(top, rate.max())
    ax.set_xlim(left=0)
    ax.set_ylim(0, top * 1.25)
    ax.set_xlabel("Dataset size (millions of residues)", fontsize=9)
    ax.set_ylabel("Throughput (10$^3$ residues/s)", fontsize=9)
    ax.tick_params(labelsize=7)
    # the GPU runs sit well above the CPU runs, so the legend goes in the gap between them
    _stream_legend(ax, sweeps, projected=False, fontsize=6.3, loc="center right")
