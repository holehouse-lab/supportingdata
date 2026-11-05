"""Subprocess worker that measures peak RAM for a single prediction run.

This module is executed as an isolated child process by
``benchmark_stream_memory.py``. Running every measurement in its own process is
deliberate: PyTorch's caching allocator and Python's own arena never release all
memory back to the OS, so measuring several dataset sizes in one long-lived
process would let a large early run inflate the peak resident set of every later
run. A fresh process per data point keeps each peak-RSS number honest.

Two modes are supported, both scoring the same synthetic FASTA file so the
numbers are directly comparable:

* ``full``   -- read the whole file into a dict and call
                :func:`metapredict.predict_disorder`, which returns every
                prediction at once (peak RAM scales with the dataset).
* ``stream`` -- iterate :func:`metapredict.predict_disorder_stream`, consuming
                and discarding each ``(header, prediction)`` tuple so nothing
                accumulates (peak RAM is bounded by ``chunk_size``).

The peak resident-set size is sampled by a background thread while the workload
runs; timing and the peak are printed to stdout as a one-line JSON record.
"""

from __future__ import annotations

import argparse
import json
import threading
import time

import psutil


def _sample_peak_rss(proc: psutil.Process, stop: threading.Event, out: dict) -> None:
    """Poll resident-set size until ``stop`` is set, recording the maximum.

    Parameters
    ----------
    proc : psutil.Process
        Handle to the current process.
    stop : threading.Event
        Sampling continues until this event is set by the main thread.
    out : dict
        Mutable container; ``out['peak']`` is updated in-place with the running
        maximum RSS in bytes.
    """
    peak = proc.memory_info().rss
    while not stop.wait(0.004):
        rss = proc.memory_info().rss
        if rss > peak:
            peak = rss
    # one final read so a spike right before shutdown is not missed
    rss = proc.memory_info().rss
    out["peak"] = max(peak, rss)


def _run_full(
    fasta_path: str, version: str, device: str, batch_size: int | None
) -> tuple[int, int]:
    """Score ``fasta_path`` in one shot with :func:`predict_disorder`.

    Returns
    -------
    tuple of int
        ``(n_sequences, n_residues)`` actually processed.
    """
    import metapredict as meta
    import protfasta

    seqs = protfasta.read_fasta(fasta_path)
    n_res = sum(len(s) for s in seqs.values())
    results = meta.predict_disorder(
        seqs, version=version, device=device, batch_size=batch_size
    )
    # Touch every result so the interpreter cannot lazily elide the work, and so
    # the full result set is genuinely resident at peak.
    total = sum(len(v[1]) for v in results.values())
    assert total == n_res
    return len(seqs), n_res


def _run_stream(
    fasta_path: str, version: str, device: str, chunk_size: int, batch_size: int | None
) -> tuple[int, int]:
    """Score ``fasta_path`` lazily with :func:`predict_disorder_stream`.

    Each yielded prediction is consumed (its length accumulated) and then dropped,
    emulating the canonical "write to disk as you go" streaming pattern so that at
    most ``chunk_size`` sequences are ever held at once.

    Returns
    -------
    tuple of int
        ``(n_sequences, n_residues)`` actually processed.
    """
    import metapredict as meta

    n_seq = 0
    n_res = 0
    for _header, (_sequence, scores) in meta.predict_disorder_stream(
        fasta_path,
        version=version,
        device=device,
        chunk_size=chunk_size,
        batch_size=batch_size,
    ):
        n_seq += 1
        n_res += len(scores)
    return n_seq, n_res


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("full", "stream"), required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--version", default="V3")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="forward-pass batch size; None uses the network/device default.",
    )
    args = parser.parse_args()

    proc = psutil.Process()

    # Warm up the network on a trivial input so that the one-off model-load cost
    # is already resident before we start sampling. Peak RSS then reflects the
    # working set of the *dataset*, not the (fixed) price of loading weights.
    import metapredict as meta

    meta.predict_disorder("MKKAADESEQ" * 3, version=args.version, device=args.device)

    baseline = proc.memory_info().rss

    sampler_out: dict = {"peak": baseline}
    stop = threading.Event()
    sampler = threading.Thread(
        target=_sample_peak_rss, args=(proc, stop, sampler_out), daemon=True
    )
    sampler.start()

    start = time.time()
    if args.mode == "full":
        n_seq, n_res = _run_full(args.fasta, args.version, args.device, args.batch_size)
    else:
        n_seq, n_res = _run_stream(
            args.fasta, args.version, args.device, args.chunk_size, args.batch_size
        )
    elapsed = time.time() - start

    stop.set()
    sampler.join()

    record = {
        "mode": args.mode,
        "version": args.version,
        "device": args.device,
        "chunk_size": args.chunk_size if args.mode == "stream" else None,
        "batch_size": args.batch_size,
        "n_sequences": n_seq,
        "n_residues": n_res,
        "elapsed_seconds": elapsed,
        "residues_per_second": n_res / elapsed if elapsed > 0 else 0.0,
        "baseline_rss_mb": baseline / 1e6,
        "peak_rss_mb": sampler_out["peak"] / 1e6,
        "delta_rss_mb": (sampler_out["peak"] - baseline) / 1e6,
    }
    print(json.dumps(record))


if __name__ == "__main__":
    main()
