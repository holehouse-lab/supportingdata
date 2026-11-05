"""Benchmark metapredict prediction time as a function of batch size (Figure S6).

For each requested network version (V1, V2, V3) and each available device (cpu,
mps, cuda), predicts disorder for the human proteome at a range of batch sizes and
records the wall-clock time. This is the measurement behind the metapredict 3.1
per-device default batch sizes: the small V1/V2 networks are dominated by
per-batch overhead on a GPU/MPS, so larger batches help a lot there.

Adapted from ``V3_manuscript/benchmark/batch_size/run_batch_size_sweep.py``. The
method is unchanged except that the batch size is now passed directly to
``predict_disorder(batch_size=...)`` (metapredict 3.1 API). The original set the
in-memory ``parameters['batch_size']``, which 3.1's per-device defaults override.

Method
------
* An off-clock warmup runs first for each (version, device), so weight loading and
  lazy MPS/CUDA initialisation are not timed.
* The device is synchronised before stopping the timer (pure inference is timed).
* Cached CUDA/MPS memory is released off the clock before each timed run, so every
  batch size starts from the same allocator state rather than inheriting blocks
  sized for the previous batch size.
* A batch size that runs out of memory (GPU, MPS or CPU) is reported on screen and
  in a summary at the end, and gets no row; the sweep carries on with the next
  batch size. Because sequences are batched longest-first, memory is set by the
  first batch, which pads every sequence to the longest one (titin, 34,350 aa).
* ``round_values=False``: for each (version, device) the smallest batch size that
  fits in memory is the reference and the maximum absolute score difference of
  every other batch size against it is recorded, confirming batch size does not
  change predictions beyond floating-point noise.

Rows are appended to a TSV in ``data/benchmarks/batch_size/`` (at the top of the
repository) as each measurement finishes.

Usage
-----
    python run_batch_size_sweep.py                          # V1-V3, all available devices
    python run_batch_size_sweep.py --devices cpu mps --machine-label "M3 Max"
    python run_batch_size_sweep.py --limit 200 --batch-sizes 32 256   # quick check
    python run_batch_size_sweep.py --cooldown 600           # 10 min cooldown per run
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

# bench_env.py (shared provenance / TSV helpers) lives one directory up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bench_env import (  # noqa: E402
    HUMAN_PROTEOME_FASTA,
    benchmark_results_dir,
    default_batch_size,
    device_available,
    host_record,
    is_out_of_memory,
    read_fasta,
    release_device_memory,
    sync_device,
    utc_now,
    utc_stamp,
)

DEFAULT_OUTPUT_DIR = benchmark_results_dir("batch_size")

# metapredict requires a power of two >= 32.
DEFAULT_BATCH_SIZES = [32, 64, 128, 256, 512, 1024]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--proteome", default=str(HUMAN_PROTEOME_FASTA), help="input FASTA file")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="directory for the results TSV")
    p.add_argument("--file-prefix", default="", help="prefix for the output filename (e.g. SMOKETEST_)")
    p.add_argument("--machine-label", default="", help="label for this machine in the figures")
    p.add_argument("--versions", nargs="+", default=["1", "2", "3"], help="network versions to test")
    p.add_argument("--devices", nargs="+", default=["cpu", "mps", "cuda"],
                   help="devices to test (unavailable ones are skipped and reported)")
    p.add_argument("--batch-sizes", nargs="+", type=int, default=DEFAULT_BATCH_SIZES)
    p.add_argument("--limit", type=int, default=None, help="use only the first N sequences (default: all)")
    p.add_argument("--warmup-n", type=int, default=200, help="sequences used for the off-clock warmup")
    p.add_argument("--cooldown", type=float, default=0.0,
                   help="seconds to sleep before each timed run so the hardware can cool "
                        "(avoids thermal throttling; e.g. 600 for a 10-minute cooldown). Default 0.")
    return p.parse_args()


def main() -> None:
    """Run the sweep and append one TSV row per (version, device, batch size)."""
    args = parse_args()
    import metapredict as meta

    sequences = read_fasta(Path(args.proteome))
    if args.limit is not None:
        sequences = dict(list(sequences.items())[: args.limit])
    warmup = dict(list(sequences.items())[: args.warmup_n])
    n_seqs = len(sequences)
    n_res = sum(len(s) for s in sequences.values())
    print(f"proteome: {Path(args.proteome).name} | {n_seqs} seqs | {n_res:,} residues")

    devices = [d for d in args.devices if device_available(d)]
    skipped = sorted(set(args.devices) - set(devices))
    if skipped:
        print(f"skipping unavailable devices: {skipped}")
    if not devices:
        raise SystemExit(f"none of the requested devices {args.devices} is available")
    print(f"devices: {devices}  |  versions: {args.versions}  |  batch sizes: {args.batch_sizes}\n")

    out_path = Path(args.output_dir) / f"{args.file_prefix}batch_size_{host_record('cpu')['hostname']}_{utc_stamp()}.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = None
    out_of_memory_runs: list[str] = []
    with open(out_path, "w", newline="") as fh:
        for version_number in args.versions:
            version = f"V{version_number}"
            for device in devices:
                host = host_record(device, args.machine_label)
                default_bs = default_batch_size(version, device)

                # off-clock warmup for this (version, device) at the default batch size
                meta.predict_disorder(warmup, version=version, device=device, show_progress_bar=False)
                sync_device(device)

                ref_time = None
                ref_scores = None
                for bs in sorted(args.batch_sizes):
                    # optional cooldown so the hardware is cool for each measurement
                    if args.cooldown:
                        print(f"  cooling down {args.cooldown:.0f}s ...", flush=True)
                        time.sleep(args.cooldown)

                    release_device_memory(device)

                    t0 = time.perf_counter()
                    try:
                        out = meta.predict_disorder(
                            sequences, version=version, device=device, batch_size=bs,
                            round_values=False, show_progress_bar=False,
                        )
                        sync_device(device)
                    except Exception as error:
                        if not is_out_of_memory(error):
                            raise
                        # no number for a batch size that does not fit; the memory
                        # the failed run held is released before the next timed run
                        out_of_memory_runs.append(f"{version} {device} bs={bs}")
                        print(f"  {version} {device:4s} bs={bs:5d}  OUT OF MEMORY - no row written", flush=True)
                        continue
                    dt = time.perf_counter() - t0

                    scores = {k: np.asarray(v[1]) for k, v in out.items()}
                    if ref_scores is None:
                        ref_time, ref_scores, max_diff = dt, scores, 0.0
                    else:
                        max_diff = max(float(np.max(np.abs(scores[k] - ref_scores[k]))) for k in scores)

                    row = {
                        "timestamp_utc": utc_now().isoformat(),
                        "network_version": version,
                        "device_used": device,
                        "batch_size": bs,
                        "default_batch_size": default_bs,
                        "n_sequences": n_seqs,
                        "n_residues": n_res,
                        "total_seconds": round(dt, 4),
                        "residues_per_second": round(n_res / dt, 1),
                        "speedup_vs_min_batch": round(ref_time / dt, 3),
                        "max_abs_diff_vs_min_batch": max_diff,
                        "cooldown_seconds": args.cooldown,
                        **host,
                    }
                    if writer is None:
                        writer = csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t")
                        writer.writeheader()
                    writer.writerow(row)
                    fh.flush()
                    print(f"  {version} {device:4s} bs={bs:5d}  {dt:8.2f}s  "
                          f"{row['residues_per_second']:>12,.0f} res/s  "
                          f"{row['speedup_vs_min_batch']:6.2f}x  diff={max_diff:.1e}", flush=True)

    print(f"\nwrote {out_path}")
    if out_of_memory_runs:
        print(f"out of memory, no row written ({len(out_of_memory_runs)}): {', '.join(out_of_memory_runs)}")


if __name__ == "__main__":
    main()
