"""Fraction of residues predicted in single-length batches (Figure S8B legend).

On x86 CPUs, PyTorch sends an LSTM batch to oneDNN's fast kernel only if every
sequence in the batch has the same length; mixed-length (packed) batches take a
slower generic path. metapredict sorts sequences by length before batching, so
the larger the dataset, the more batches hold a single length. This is why
``predict_disorder`` speeds up with dataset size on the Xeon CPU in Figure S8B
while ``predict_disorder_stream``, which sorts each chunk on its own, does not.

This script reproduces the fractions quoted in the Figure S8 legend. It needs no
prediction and no FASTA files: it redraws the sequence lengths of each synthetic
dataset from the benchmark's seed and counts the residues that fall in batches
whose sequences all share one length.

Usage
-----
    python single_length_batch_fraction.py
"""

from __future__ import annotations

import numpy as np

from benchmark_stream_memory import (
    FIXED_BATCH_SIZE,
    FIXED_CHUNK_SIZE,
    FULL_SIZES,
    RNG_SEED,
    draw_lengths,
)


def single_length_residues(lengths: np.ndarray, batch_size: int) -> int:
    """Count the residues that fall in batches made of one sequence length.

    Sequences are sorted longest-first and cut into consecutive batches of
    ``batch_size``, as metapredict does within one prediction call.

    Parameters
    ----------
    lengths : numpy.ndarray
        Sequence lengths, in residues, of the sequences predicted in one call.
    batch_size : int
        Number of sequences per batch.

    Returns
    -------
    int
        Number of residues in batches whose sequences all have the same length.
    """
    assert lengths.size > 0, "no sequences to batch"
    ordered = np.sort(lengths)[::-1]
    residues = 0
    for start in range(0, ordered.size, batch_size):
        batch = ordered[start : start + batch_size]
        # sorted, so the batch has a single length if its ends are equal
        if batch[0] == batch[-1]:
            residues += int(batch.sum())
    return residues


def single_length_fractions(n_sequences: int) -> tuple[int, float, float]:
    """Single-length-batch fractions for one synthetic dataset.

    Parameters
    ----------
    n_sequences : int
        Number of sequences in the dataset (one of the benchmark's sweep sizes).

    Returns
    -------
    tuple
        ``(n_residues, full_fraction, stream_fraction)``: the dataset size in
        residues, and the fraction of residues in single-length batches when the
        whole dataset is sorted at once (``predict_disorder``) and when each
        chunk of :data:`FIXED_CHUNK_SIZE` sequences is sorted on its own
        (``predict_disorder_stream``).
    """
    lengths = draw_lengths(np.random.default_rng(RNG_SEED), n_sequences)
    n_residues = int(lengths.sum())
    full_residues = single_length_residues(lengths, FIXED_BATCH_SIZE)
    stream_residues = 0
    for start in range(0, n_sequences, FIXED_CHUNK_SIZE):
        chunk = lengths[start : start + FIXED_CHUNK_SIZE]
        stream_residues += single_length_residues(chunk, FIXED_BATCH_SIZE)
    return n_residues, full_residues / n_residues, stream_residues / n_residues


def main() -> None:
    """Print the fractions for every dataset of the full stream sweep."""
    print(f"batch size {FIXED_BATCH_SIZE}, chunk size {FIXED_CHUNK_SIZE:,}")
    print("sequences\tresidues\tpredict_disorder\tpredict_disorder_stream")
    for n_sequences in FULL_SIZES:
        n_residues, full_fraction, stream_fraction = single_length_fractions(
            n_sequences
        )
        print(
            f"{n_sequences:,}\t{n_residues:,}\t{full_fraction:.1%}\t{stream_fraction:.1%}"
        )


if __name__ == "__main__":
    main()
