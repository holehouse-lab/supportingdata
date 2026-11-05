"""Time batch vs one-at-a-time (sequential) disorder prediction (Figure S7).

Draws random subsets of increasing size from the human proteome and times
predicting them either one sequence at a time (a loop of single-string
``predict_disorder`` calls, the way a naive script would) or all at once (a single
batched call on the list), for each metapredict network.

Converted from the 2025 notebook
``supplemental_figure_NEW_batch_vs_individual.ipynb``: same subset sizes and the
same two timing functions, now looped over V1/V2/V3, seeded, and written to TSV.

Notes
-----
* metapredict always predicts a single string on the CPU, so sequential rows are
  CPU timings whatever ``--device`` is; batch rows use ``--device`` (default:
  metapredict's own per-network choice).
* Batched calls use metapredict's default batch size for that network/device.
* Each subset size is drawn once per replicate with a seeded RNG, and the same
  subset is timed both ways.

Usage
-----
    python run_batch_vs_sequential.py                      # V1-V3, default device
    python run_batch_vs_sequential.py --device cpu --machine-label "x86-64"
    python run_batch_vs_sequential.py --counts 10 50 100 --versions 3   # quick check
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

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
    sync_device,
    utc_now,
    utc_stamp,
    write_tsv,
)

DEFAULT_OUTPUT_DIR = benchmark_results_dir("batch_vs_sequential")

# Subset sizes (number of sequences) from the original notebook.
DEFAULT_COUNTS = [10, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 8000, 10000, 15000]
DEFAULT_SEED = 20250614  # date of the original notebook


def sequential(seqs: list[str], version: str) -> float:
    """Predict each sequence with its own call; return elapsed seconds.

    Parameters
    ----------
    seqs : list of str
        Sequences to predict.
    version : str
        Network version (``"V1"``, ``"V2"``, ``"V3"``).

    Returns
    -------
    float
        Wall-clock seconds for the whole loop.
    """
    import metapredict as meta

    start = time.perf_counter()
    for s in seqs:
        meta.predict_disorder(s, version=version)
    return time.perf_counter() - start


def batch(seqs: list[str], version: str, device: str) -> float:
    """Predict all sequences in one batched call; return elapsed seconds.

    Parameters
    ----------
    seqs : list of str
        Sequences to predict.
    version : str
        Network version (``"V1"``, ``"V2"``, ``"V3"``).
    device : str
        torch device for the batched call.

    Returns
    -------
    float
        Wall-clock seconds, with the device synchronised before stopping the clock.
    """
    import metapredict as meta

    start = time.perf_counter()
    meta.predict_disorder(seqs, version=version, device=device, show_progress_bar=False)
    sync_device(device)
    return time.perf_counter() - start


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--proteome", default=str(HUMAN_PROTEOME_FASTA), help="input FASTA file")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="directory for the results TSV")
    p.add_argument("--file-prefix", default="", help="prefix for the output filename (e.g. SMOKETEST_)")
    p.add_argument("--machine-label", default="", help="label for this machine in the figures")
    p.add_argument("--versions", nargs="+", default=["1", "2", "3"], help="network versions to test")
    p.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto",
                   help="device for the batched calls ('auto' = metapredict's per-network default)")
    p.add_argument("--counts", nargs="+", type=int, default=DEFAULT_COUNTS, help="subset sizes")
    p.add_argument("--replicates", type=int, default=1, help="independent random subsets per size")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED, help="RNG seed for subset selection")
    return p.parse_args()


def main() -> None:
    """Run the sweep and write one TSV row per (version, mode, size, replicate)."""
    args = parse_args()
    proteome = read_fasta(Path(args.proteome))
    headers = list(proteome.keys())
    if max(args.counts) > len(headers):
        raise SystemExit(f"cannot sample {max(args.counts)} sequences from {len(headers)}")
    rng = random.Random(args.seed)
    out_path = None
    hosts: dict[str, dict] = {}  # provenance per device, looked up once

    rows = []
    for version_number in args.versions:
        version = f"V{version_number}"
        device = default_device(version) if args.device == "auto" else args.device
        if not device_available(device):
            raise SystemExit(f"device {device} is not available on this machine")

        # off-clock warmup so weight loading / device init is not timed
        warm = [proteome[h] for h in headers[:20]]
        sequential(warm[:2], version)
        batch(warm, version, device)

        for replicate in range(args.replicates):
            for count in args.counts:
                seqs = [proteome[h] for h in rng.sample(headers, count)]
                n_res = sum(len(s) for s in seqs)
                timings = {"sequential": (sequential(seqs, version), "cpu", ""),
                           "batch": (batch(seqs, version, device), device,
                                     default_batch_size(version, device))}
                for mode, (seconds, mode_device, batch_size) in timings.items():
                    if mode_device not in hosts:
                        hosts[mode_device] = host_record(mode_device, args.machine_label)
                    rows.append({
                        "timestamp_utc": utc_now().isoformat(),
                        "network_version": version,
                        "mode": mode,
                        "device_used": mode_device,
                        "batch_size": batch_size,
                        "n_sequences": count,
                        "n_residues": n_res,
                        "replicate": replicate,
                        "seed": args.seed,
                        "total_seconds": round(seconds, 4),
                        **hosts[mode_device],
                    })
                print(f"  {version} n={count:6d}  sequential {timings['sequential'][0]:8.2f}s  "
                      f"batch ({device}) {timings['batch'][0]:7.2f}s", flush=True)

        # rewrite the TSV after each version so an interrupted run keeps its data
        if out_path is None:
            out_path = (Path(args.output_dir)
                        / f"{args.file_prefix}batch_vs_sequential_{rows[0]['hostname']}_{utc_stamp()}.tsv")
        write_tsv(out_path, rows)


if __name__ == "__main__":
    main()
