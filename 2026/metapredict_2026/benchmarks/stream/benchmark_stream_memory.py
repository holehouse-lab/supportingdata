"""Benchmark the out-of-core advantage of predict_disorder_stream over predict_disorder.

This is the Figure S8 benchmark (originally ``figures/stream_performance``). The
experiment is unchanged; this copy writes a TSV in ``data/benchmarks/stream/`` (one row per
measurement, with metapredict and hardware provenance) instead of JSON.

Motivation
----------
:func:`metapredict.predict_disorder` loads a whole dataset, predicts it, and
returns *every* result at once, so the entire input dictionary and the entire
output dictionary are co-resident in RAM for the duration of the call. Its peak
memory therefore grows with the size of the dataset, and a file whose sequences
(plus their per-residue scores) do not fit in RAM cannot be scored at all.

:func:`metapredict.predict_disorder_stream` reads the FASTA lazily and predicts a
fixed number of sequences at a time (``chunk_size``), yielding each result and
then dropping it. At most ``chunk_size`` sequences are ever held, so its peak
memory is *flat* in the dataset size and it can score arbitrarily large files.

Peak RSS is dominated by the forward-pass activation memory of the most-padded
batch, which is set by ``batch_size`` (``batch_size x max_seq_len x hidden``),
**not** by ``chunk_size``. To expose the streaming-vs-full difference -- the
co-resident input/output data -- ``batch_size`` is held fixed and small for *both*
functions so the activation floor is low and constant, and only the dataset size
varies.

The experiment is fully synthetic (fixed RNG seed, no downloads). Every
measurement runs in a fresh subprocess (``stream_memory_worker.py``) so peak-RSS
numbers cannot leak between runs. Rows are appended to
``data/benchmarks/stream/stream_<host>_<version>_<device>_bs<batch>_chunk<chunk>[_quick].tsv``; a
re-run with the same configuration resumes from that file, skipping points that
are already measured.

Usage
-----
    python benchmark_stream_memory.py --quick   # 3 sizes, a few minutes
    python benchmark_stream_memory.py           # 7 sizes up to ~1e9 residues, ~4 h on an M3 Max
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

# bench_env.py (shared provenance / TSV helpers) lives one directory up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bench_env import BENCHMARK_DATA_DIR, benchmark_results_dir, host_record, utc_now  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = str(benchmark_results_dir("stream"))
# Synthetic FASTA inputs (~2 GB for the full sweep; rebuilt from a fixed seed, not
# tracked). Set METAPREDICT_STREAM_DATA_DIR to write them elsewhere, e.g. a local
# scratch disk outside a synced folder.
DATA_DIR = os.environ.get("METAPREDICT_STREAM_DATA_DIR",
                          str(BENCHMARK_DATA_DIR / "stream_synthetic_data"))
WORKER = os.path.join(HERE, "stream_memory_worker.py")

# Reproducible synthetic-proteome parameters. Sequence lengths are drawn from a
# log-normal distribution loosely matched to real proteomes (median ~300 aa) and
# clamped to a sane range; residues are uniform over the 20 standard amino acids.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
LENGTH_LOGNORM_MEAN = 5.7  # exp(5.7) ~= 300 residues
LENGTH_LOGNORM_SIGMA = 0.55
LENGTH_MIN = 40
LENGTH_MAX = 2500
RNG_SEED = 20260705

# Fixed low batch size for BOTH functions: keeps the (constant) activation floor
# small so the dataset-scaling term is what varies. See module docstring.
FIXED_BATCH_SIZE = 32
# Fixed streaming chunk: the number of sequences held at once by the streaming
# call. Its exact value barely matters for peak RAM (that is the whole point).
# This is metapredict 3.1's default; the original sweep used 5,000.
FIXED_CHUNK_SIZE = 20000

# Sweep sizes (number of sequences). The full sweep spans ~7M -> ~1000M residues,
# i.e. from below the human proteome (~11.5M) to ~85x it; the top point is ~1 h of
# CPU compute per function.
QUICK_SIZES = [20000, 80000, 200000]
FULL_SIZES = [20000, 80000, 250000, 500000, 720000, 1440000, 2880000]

# worker memory is in decimal MB (1e6 bytes); results are stored in GiB (2**30 bytes)
MB_TO_GIB = 1e6 / 2**30


def make_fasta(n_sequences: int, path: str) -> tuple[int, int]:
    """Write a reproducible synthetic FASTA file of ``n_sequences`` sequences.

    The RNG is re-seeded from :data:`RNG_SEED` on every call, so a given
    ``n_sequences`` always produces byte-identical output and cached files are
    reused across re-runs.

    Parameters
    ----------
    n_sequences : int
        Number of sequences to generate.
    path : str
        Destination FASTA path.

    Returns
    -------
    tuple of int
        ``(n_sequences, n_residues)`` written.
    """
    rng = np.random.default_rng(RNG_SEED)
    lengths = np.clip(
        rng.lognormal(LENGTH_LOGNORM_MEAN, LENGTH_LOGNORM_SIGMA, n_sequences).astype(
            int
        ),
        LENGTH_MIN,
        LENGTH_MAX,
    )
    alphabet = np.frombuffer(AA_ALPHABET.encode(), dtype="S1")

    # Vectorised generation: draw every residue for the whole file in one shot
    # (as int8 indices, ~1 byte/residue, so even a billion residues is ~1 GB of
    # transient), then slice per sequence and write bytes directly. This is orders
    # of magnitude faster than a per-sequence rng.choice loop at multi-million
    # sequence counts.
    total = int(lengths.sum())
    idx = rng.integers(0, len(alphabet), size=total, dtype=np.int8)
    residues = alphabet[idx]  # S1 array of length `total`
    ends = np.cumsum(lengths)
    starts = ends - lengths

    with open(path, "wb") as handle:
        buf = bytearray()
        for i in range(n_sequences):
            handle_seq = residues[starts[i] : ends[i]].tobytes()
            buf += b">synthetic_%d\n%s\n" % (i, handle_seq)
            if len(buf) >= 8 << 20:  # flush every ~8 MB
                handle.write(buf)
                buf = bytearray()
        if buf:
            handle.write(buf)
    return n_sequences, total


def ensure_fasta(n_sequences: int) -> str:
    """Return the path to the synthetic FASTA for ``n_sequences``, building it once."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, f"synthetic_{n_sequences}.fasta")
    if not os.path.exists(path):
        _, n_res = make_fasta(n_sequences, path)
        print(f"  built {path} ({n_sequences:,} seqs, {n_res:,} residues)")
    return path


def run_worker(
    mode: str,
    fasta: str,
    version: str,
    device: str,
    batch_size: int,
    chunk_size: int = FIXED_CHUNK_SIZE,
) -> dict[str, Any]:
    """Invoke ``stream_memory_worker.py`` in a subprocess and parse its JSON record."""
    cmd = [
        sys.executable,
        WORKER,
        "--mode",
        mode,
        "--fasta",
        fasta,
        "--version",
        version,
        "--device",
        device,
        "--batch-size",
        str(batch_size),
        "--chunk-size",
        str(chunk_size),
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"worker failed ({mode}, n from {fasta}):\n{completed.stderr}"
        )
    line = [ln for ln in completed.stdout.splitlines() if ln.strip()][-1]
    return json.loads(line)




def results_path(args: argparse.Namespace, hostname: str, sizes_tag: str) -> str:
    """TSV path for one sweep configuration (re-runs with the same config resume it)."""
    return os.path.join(
        args.output_dir,
        f"{args.file_prefix}stream_{hostname}_{args.version}_{args.device}_"
        f"bs{args.batch_size}_chunk{FIXED_CHUNK_SIZE}{sizes_tag}.tsv",
    )


def main() -> None:
    """Run the dataset-size sweep for both functions, appending TSV rows."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="V3", help="metapredict network version")
    parser.add_argument(
        "--device",
        default="cpu",
        help="torch device. 'cpu' (default) gives a clean host-RAM signal and, "
        f"at batch_size={FIXED_BATCH_SIZE}, is actually faster than 'mps' here "
        "(a 32-wide batch is too small to amortise GPU dispatch overhead).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=FIXED_BATCH_SIZE,
        help="forward-pass batch size, held fixed for both functions.",
    )
    parser.add_argument("--quick", action="store_true", help="small, fast sweep (3 sizes).")
    parser.add_argument("--sizes", nargs="+", type=int, default=None,
                        help="explicit sequence counts (overrides --quick; for smoke tests).")
    parser.add_argument("--machine-label", default="", help="label for this machine in the figures")
    parser.add_argument("--output-dir", default=OUTPUT_DIR, help="directory for the results TSV")
    parser.add_argument("--file-prefix", default="", help="prefix for the output filename")
    args = parser.parse_args()

    if args.sizes:
        scaling_sizes, sizes_tag = args.sizes, "_custom"
    elif args.quick:
        scaling_sizes, sizes_tag = QUICK_SIZES, "_quick"
    else:
        scaling_sizes, sizes_tag = FULL_SIZES, ""

    os.makedirs(args.output_dir, exist_ok=True)
    host = host_record(args.device, args.machine_label)
    out = results_path(args, host["hostname"], sizes_tag)

    # Resume: reuse already-measured points from this configuration's TSV, so a
    # long sweep can be interrupted and restarted without repeating work.
    done: set[tuple[str, int]] = set()
    if os.path.exists(out):
        with open(out) as handle:
            done = {(r["mode"], int(r["n_sequences"])) for r in csv.DictReader(handle, delimiter="\t")}
        print(f"resuming: {len(done)} measurement(s) already in {out}")

    print(
        f"device = {args.device}, version = {args.version}, "
        f"batch_size = {args.batch_size}, chunk_size = {FIXED_CHUNK_SIZE}"
    )
    print("\nscaling sweep: peak RAM vs dataset size (predict_disorder vs stream)")

    for n in scaling_sizes:
        fasta = None
        for mode in ("full", "stream"):
            if (mode, n) in done:
                print(f"  {mode:6s} n={n:>9,}  (cached, skipping)")
                continue
            if fasta is None:
                fasta = ensure_fasta(n)
            rec = run_worker(mode, fasta, args.version, args.device, args.batch_size)
            row: dict[str, Any] = {
                "timestamp_utc": utc_now().isoformat(),
                "mode": rec["mode"],
                "network_version": rec["version"],
                "device_used": rec["device"],
                "batch_size": rec["batch_size"],
                "chunk_size": FIXED_CHUNK_SIZE,
                "sweep": sizes_tag.strip("_") or "full",
                "rng_seed": RNG_SEED,
                "n_sequences": rec["n_sequences"],
                "n_residues": rec["n_residues"],
                "total_seconds": round(rec["elapsed_seconds"], 4),
                "residues_per_second": round(rec["residues_per_second"], 1),
                # the worker reports decimal MB (bytes / 1e6); store GiB like the
                # proteome benchmark
                "baseline_rss_gb": round(rec["baseline_rss_mb"] * MB_TO_GIB, 4),
                "peak_rss_gb": round(rec["peak_rss_mb"] * MB_TO_GIB, 4),
                "delta_rss_gb": round(rec["delta_rss_mb"] * MB_TO_GIB, 4),
                **host,
            }
            # append, writing the header only for a new file (persist every point)
            new_file = not os.path.exists(out)
            with open(out, "a", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
                if new_file:
                    writer.writeheader()
                writer.writerow(row)
            print(
                f"  {mode:6s} n={rec['n_sequences']:>9,} "
                f"res={rec['n_residues']:>12,} "
                f"peak={rec['peak_rss_mb']:8.0f} MB "
                f"{rec['residues_per_second'] / 1e3:8.0f}k res/s "
                f"({rec['elapsed_seconds']:.0f} s)"
            )

    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
