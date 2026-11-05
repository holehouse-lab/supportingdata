"""Compute and memory benchmark for metapredict disorder prediction over a proteome.

Benchmarks one metapredict network (V1, V2 or V3) on one device (CPU, Apple MPS
or CUDA GPU) while predicting disorder for every sequence in a FASTA file (by
default the human reference proteome, UP000005640), and writes one TSV row
describing the run into ``data/benchmarks/proteome/`` (at the top of the repository).

This is the Figure 1E / Figure S5 benchmark (originally
``figures/Figure_1E/benchmark_metapredict.py``). The resource monitor and timed
region are unchanged; this copy writes a flat TSV row instead of JSON, records the
exact metapredict install (version, path, git commit) and the batch size it
resolved, and adds ``--device auto`` to time metapredict's own device choice.

* File I/O and the one-time model / device initialisation (an untimed warmup)
  are excluded from the timed region; GPU work is synchronised before stopping
  the clock.
* ``batch_size`` is deliberately NOT passed to metapredict, so each run uses the
  metapredict 3.1 per-network, per-device default (256 on CPU/CUDA, 512 on MPS).
* Peak RSS, CPU utilisation and device memory are sampled by a background thread
  during the timed region.

Examples
--------
Benchmark V3 on the Mac GPU::

    python benchmark_metapredict.py --version 3 --device mps

Benchmark V1 on a Linux CPU, labelling the machine for the figures::

    python benchmark_metapredict.py --version 1 --device cpu --machine-label "x86-64"
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import psutil

# bench_env.py (shared provenance / TSV helpers) lives one directory up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bench_env import (  # noqa: E402
    HUMAN_PROTEOME_FASTA,
    benchmark_results_dir,
    default_batch_size,
    default_device,
    device_available,
    host_record,
    read_fasta,
    utc_now,
    utc_stamp,
    write_tsv,
)

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

DEFAULT_OUTPUT_DIR = benchmark_results_dir("proteome")

# Map the user-facing device names onto the strings metapredict / torch expect.
DEVICE_MAP: dict[str, str] = {"cpu": "cpu", "mps": "mps", "gpu": "cuda"}

# How often (seconds) the resource monitor samples RSS / CPU / GPU memory.
MONITOR_INTERVAL_S = 0.1

# Number of sequences used for the untimed warmup pass (loads weights on device).
WARMUP_N_SEQUENCES = 4

MB_PER_GB = 1024.0




# --------------------------------------------------------------------------- #
# Resource monitoring
# --------------------------------------------------------------------------- #


@dataclass
class MonitorSummary:
    """Container for the aggregated output of a :class:`ResourceMonitor` run.

    Attributes
    ----------
    baseline_rss_mb : float
        Process resident set size, in mebibytes, sampled immediately before the
        timed region began.
    peak_rss_mb : float
        Maximum process resident set size, in mebibytes, observed during the
        timed region.
    mean_rss_mb : float
        Mean process resident set size, in mebibytes, across all samples.
    peak_gpu_mb : Optional[float]
        Peak device-side memory, in mebibytes. For CUDA this is
        ``torch.cuda.max_memory_allocated``; for MPS it is the maximum sampled
        ``torch.mps.driver_allocated_memory``. ``None`` for CPU-only runs.
    cpu_percent_mean : float
        Mean process CPU utilisation across samples. Values may exceed 100 as
        they aggregate over all logical cores used by the process.
    cpu_percent_max : float
        Maximum process CPU utilisation observed.
    n_samples : int
        Number of samples collected.
    """

    baseline_rss_mb: float
    peak_rss_mb: float
    mean_rss_mb: float
    peak_gpu_mb: Optional[float]
    cpu_percent_mean: float
    cpu_percent_max: float
    n_samples: int


class ResourceMonitor:
    """Background sampler of process memory, CPU, and device memory.

    The monitor runs in a daemon thread and records the process resident set
    size (RSS), the process CPU utilisation, and — when a GPU device is in use —
    the device-side allocated memory, at a fixed interval. Aggregated statistics
    are retrieved via :meth:`summarise` after :meth:`stop`.

    Parameters
    ----------
    device : str
        The torch device string in use (``"cpu"``, ``"mps"`` or ``"cuda"``).
        Governs whether and how GPU memory is sampled.
    interval_s : float, optional
        Sampling interval in seconds. Default is :data:`MONITOR_INTERVAL_S`.
    """

    def __init__(self, device: str, interval_s: float = MONITOR_INTERVAL_S) -> None:
        self._device = device
        self._interval_s = interval_s
        self._proc = psutil.Process(os.getpid())
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

        self._rss_samples: list[float] = []
        self._cpu_samples: list[float] = []
        self._gpu_samples: list[float] = []
        self._baseline_rss_mb: float = self._rss_mb()

        # torch is only needed for GPU memory sampling; import lazily.
        self._torch = None
        if self._device in ("cuda", "mps"):
            import torch  # local import keeps CPU runs torch-free at monitor level

            self._torch = torch

    # -- sampling helpers --------------------------------------------------- #

    def _rss_mb(self) -> float:
        """Return the current process RSS in mebibytes."""
        return self._proc.memory_info().rss / (1024.0 * 1024.0)

    def _gpu_mb(self) -> Optional[float]:
        """Return current device-allocated memory in mebibytes, or ``None``."""
        if self._torch is None:
            return None
        try:
            if self._device == "mps":
                return self._torch.mps.driver_allocated_memory() / (1024.0 * 1024.0)
            if self._device == "cuda":
                return self._torch.cuda.memory_allocated() / (1024.0 * 1024.0)
        except Exception:
            return None
        return None

    # -- thread lifecycle --------------------------------------------------- #

    def _run(self) -> None:
        # Prime CPU percent measurement; the first call always returns 0.0.
        self._proc.cpu_percent(interval=None)
        while not self._stop_event.is_set():
            self._rss_samples.append(self._rss_mb())
            self._cpu_samples.append(self._proc.cpu_percent(interval=None))
            gpu = self._gpu_mb()
            if gpu is not None:
                self._gpu_samples.append(gpu)
            self._stop_event.wait(self._interval_s)

    def start(self) -> None:
        """Record the memory baseline and start the sampling thread."""
        self._baseline_rss_mb = self._rss_mb()
        if self._device == "cuda" and self._torch is not None:
            # Reset the CUDA peak counter so max_memory_allocated reflects only
            # the timed region.
            self._torch.cuda.reset_peak_memory_stats()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        """Signal the sampling thread to stop and join it."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join()

    def summarise(self) -> MonitorSummary:
        """Aggregate the collected samples into a :class:`MonitorSummary`.

        Returns
        -------
        MonitorSummary
            The aggregated memory / CPU / GPU statistics for the run.
        """
        rss = self._rss_samples or [self._baseline_rss_mb]
        cpu = self._cpu_samples or [0.0]

        peak_gpu_mb: Optional[float] = None
        if self._device == "cuda" and self._torch is not None:
            try:
                peak_gpu_mb = self._torch.cuda.max_memory_allocated() / (1024.0 * 1024.0)
            except Exception:
                peak_gpu_mb = max(self._gpu_samples) if self._gpu_samples else None
        elif self._device == "mps" and self._gpu_samples:
            peak_gpu_mb = max(self._gpu_samples)

        return MonitorSummary(
            baseline_rss_mb=round(self._baseline_rss_mb, 3),
            peak_rss_mb=round(max(rss), 3),
            mean_rss_mb=round(sum(rss) / len(rss), 3),
            peak_gpu_mb=round(peak_gpu_mb, 3) if peak_gpu_mb is not None else None,
            cpu_percent_mean=round(sum(cpu) / len(cpu), 2),
            cpu_percent_max=round(max(cpu), 2),
            n_samples=len(rss),
        )


# --------------------------------------------------------------------------- #
# The benchmark itself
# --------------------------------------------------------------------------- #


@dataclass
class BenchmarkResult:
    """Timing and throughput results for a single benchmark run.

    Attributes
    ----------
    total_seconds : float
        Wall-clock seconds spent inside the timed inference region (excludes
        file I/O and warmup / model loading).
    residues_per_second : float
        Total residues predicted divided by ``total_seconds``.
    sequences_per_second : float
        Total sequences predicted divided by ``total_seconds``.
    warmup_seconds : float
        Wall-clock seconds spent on the untimed warmup pass (model / device
        initialisation). Reported for transparency but excluded from throughput.
    n_sequences : int
        Number of sequences predicted.
    n_residues : int
        Total number of residues predicted.
    """

    total_seconds: float
    residues_per_second: float
    sequences_per_second: float
    warmup_seconds: float
    n_sequences: int
    n_residues: int


def _sync_device(device: str) -> None:
    """Block until all queued work on the given device has completed.

    Ensures GPU timing captures the full compute cost rather than only the time
    to enqueue asynchronous kernels.

    Parameters
    ----------
    device : str
        The torch device string (``"cpu"``, ``"mps"`` or ``"cuda"``).
    """
    if device == "cuda":
        import torch

        if torch.cuda.is_available():
            torch.cuda.synchronize()
    elif device == "mps":
        import torch

        if torch.backends.mps.is_available():
            torch.mps.synchronize()


def run_benchmark(
    sequences: dict[str, str],
    *,
    version: str,
    device: str,
    domain_decomposition: bool,
    do_warmup: bool = True,
) -> tuple[BenchmarkResult, MonitorSummary]:
    """Run the timed metapredict prediction over all sequences.

    Parameters
    ----------
    sequences : dict of str to str
        Sequences to predict, as a header -> sequence mapping.
    version : str
        metapredict network version, one of ``"V1"``, ``"V2"``, ``"V3"``.
    device : str
        torch device string: ``"cpu"``, ``"mps"`` or ``"cuda"``.
    domain_decomposition : bool
        If ``True``, request full IDR / folded-domain decomposition
        (``return_domains=True``); otherwise only raw per-residue disorder
        scores are computed.
    do_warmup : bool, optional
        If ``True`` (default), run an untimed warmup pass on a small batch to
        load model weights onto the device before timing.

    Returns
    -------
    BenchmarkResult
        Timing and throughput metrics for the run.
    MonitorSummary
        Memory and CPU statistics gathered during the timed region.
    """
    from metapredict import meta

    n_sequences = len(sequences)
    n_residues = sum(len(s) for s in sequences.values())

    # --- Untimed warmup: force weight load / lazy device init off the clock. --
    warmup_seconds = 0.0
    if do_warmup:
        warmup_keys = list(sequences.keys())[:WARMUP_N_SEQUENCES]
        warmup_batch = {k: sequences[k] for k in warmup_keys}
        w0 = time.perf_counter()
        meta.predict_disorder(
            warmup_batch,
            version=version,
            device=device,
            return_domains=domain_decomposition,
            show_progress_bar=False,
        )
        _sync_device(device)
        warmup_seconds = time.perf_counter() - w0

    # --- Timed region. ----------------------------------------------------- #
    monitor = ResourceMonitor(device=device)
    monitor.start()
    t0 = time.perf_counter()

    meta.predict_disorder(
        sequences,
        version=version,
        device=device,
        return_domains=domain_decomposition,
        show_progress_bar=False,
    )

    _sync_device(device)
    total_seconds = time.perf_counter() - t0
    monitor.stop()

    monitor_summary = monitor.summarise()

    result = BenchmarkResult(
        total_seconds=round(total_seconds, 4),
        residues_per_second=round(n_residues / total_seconds, 2),
        sequences_per_second=round(n_sequences / total_seconds, 2),
        warmup_seconds=round(warmup_seconds, 4),
        n_sequences=n_sequences,
        n_residues=n_residues,
    )
    return result, monitor_summary




# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


def resolve_device(user_device: str, version: str) -> str:
    """Translate and validate a user device choice into a torch device string.

    Parameters
    ----------
    user_device : str
        One of ``"cpu"``, ``"mps"``, ``"gpu"`` or ``"auto"`` (metapredict's own
        per-network choice when no device is given).
    version : str
        Network version (``"V1"``/``"V2"``/``"V3"``), used for ``"auto"``.

    Returns
    -------
    str
        The torch device string (``"cpu"``, ``"mps"``, ``"cuda"``).

    Raises
    ------
    SystemExit
        If the requested accelerator is not available on this machine, so the
        benchmark fails loudly rather than silently falling back to CPU.
    """
    if user_device == "auto":
        return default_device(version)
    device = DEVICE_MAP[user_device]
    if not device_available(device):
        raise SystemExit(f"Requested device '{user_device}' ({device}) is not available on this machine.")
    return device


def build_row(
    *,
    args: argparse.Namespace,
    fasta_path: Path,
    device_used: str,
    result: BenchmarkResult,
    memory: MonitorSummary,
) -> dict:
    """Assemble the flat TSV record for one benchmark run.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command-line arguments.
    fasta_path : pathlib.Path
        FASTA file that was benchmarked.
    device_used : str
        torch device string actually used.
    result : BenchmarkResult
        Timing / throughput results.
    memory : MonitorSummary
        Memory / CPU results.

    Returns
    -------
    dict
        One row: run metadata, hardware and metapredict provenance, timing and
        memory (memory in GB).
    """
    version = f"V{args.version}"
    row = {
        "timestamp_utc": utc_now().isoformat(),
        "network_version": version,
        "device_requested": args.device,
        "device_used": device_used,
        "batch_size": default_batch_size(version, device_used),
        "domain_decomposition": bool(args.domain_decomposition),
        "fasta": fasta_path.name,
        "sequence_limit": args.limit if args.limit is not None else "",
        "n_sequences": result.n_sequences,
        "n_residues": result.n_residues,
        "warmup_performed": not args.no_warmup,
        "warmup_seconds": result.warmup_seconds,
        "total_seconds": result.total_seconds,
        "residues_per_second": result.residues_per_second,
        "sequences_per_second": result.sequences_per_second,
        "baseline_rss_gb": round(memory.baseline_rss_mb / MB_PER_GB, 4),
        "peak_rss_gb": round(memory.peak_rss_mb / MB_PER_GB, 4),
        "mean_rss_gb": round(memory.mean_rss_mb / MB_PER_GB, 4),
        "peak_gpu_gb": round(memory.peak_gpu_mb / MB_PER_GB, 4) if memory.peak_gpu_mb is not None else "",
        "cpu_percent_mean": memory.cpu_percent_mean,
        "cpu_percent_max": memory.cpu_percent_max,
        "monitor_samples": memory.n_samples,
    }
    row.update(host_record(device_used, args.machine_label))
    return row


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark metapredict (V1/V2/V3) disorder prediction over a proteome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", type=int, choices=(1, 2, 3), required=True,
                        help="metapredict network version to benchmark.")
    parser.add_argument("--device", choices=("cpu", "mps", "gpu", "auto"), required=True,
                        help="Compute device. 'gpu' maps to CUDA; 'auto' uses metapredict's default.")
    parser.add_argument("--machine-label", default="",
                        help="Label for this machine in the figures (e.g. 'x86-64'); "
                             "defaults to the CPU/GPU name.")
    parser.add_argument("--domain-decomposition", action="store_true",
                        help="Return full IDR / folded-domain decomposition (return_domains=True).")
    parser.add_argument("--fasta", default=str(HUMAN_PROTEOME_FASTA), help="Input FASTA file.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR),
                        help="Directory for the results TSV.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Benchmark only the first N sequences (quick tests).")
    parser.add_argument("--no-warmup", action="store_true",
                        help="Skip the untimed warmup pass (weight loading will be timed).")
    parser.add_argument("--file-prefix", default="",
                        help="Prefix for the output filename (e.g. SMOKETEST_).")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    """Run one benchmark and write a one-row TSV into the output directory."""
    args = parse_args(argv)

    fasta_path = Path(args.fasta)
    if not fasta_path.is_file():
        raise SystemExit(f"FASTA file not found: {fasta_path}")

    version = f"V{args.version}"
    device_used = resolve_device(args.device, version)

    print(f"[benchmark] metapredict {version} on device '{device_used}' "
          f"(default batch size {default_batch_size(version, device_used)})")
    sequences = read_fasta(fasta_path)
    if args.limit is not None:
        keys = list(sequences.keys())[: args.limit]
        sequences = {k: sequences[k] for k in keys}
    print(f"[benchmark] {len(sequences)} sequences, {sum(len(s) for s in sequences.values())} residues")

    result, memory = run_benchmark(
        sequences,
        version=version,
        device=device_used,
        domain_decomposition=args.domain_decomposition,
        do_warmup=not args.no_warmup,
    )
    row = build_row(args=args, fasta_path=fasta_path, device_used=device_used,
                    result=result, memory=memory)

    domain_tag = "domains" if args.domain_decomposition else "scores"
    name = (f"{args.file_prefix}proteome_{version}_{device_used}_{domain_tag}_"
            f"{row['hostname']}_{utc_stamp()}.tsv")
    write_tsv(Path(args.output_dir) / name, [row])

    print(f"[benchmark] {result.total_seconds} s, {result.residues_per_second:,.0f} residues/s, "
          f"peak RSS {row['peak_rss_gb']} GB, peak GPU {row['peak_gpu_gb'] or 'n/a'} GB")


if __name__ == "__main__":
    main()
