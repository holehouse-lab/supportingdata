"""Load the human proteome / transcription factor IDR tables for Figure S16.

Replaces the SHEPHARD Proteome objects used in the original notebook with plain
dictionaries: sequences come from the UniProt FASTA files and IDRs from the domain
tables written by ``build_idr_tables.py``.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np

from build_idr_tables import (HUMAN_PROTEOME_FASTA, HUMAN_PROTEOME_IDRS, HUMAN_TF_FASTA,
                              HUMAN_TF_IDRS, read_uniprot_fasta)

AMINO_ACIDS = ['A', 'L', 'V', 'I', 'M', 'F', 'Y', 'W', 'P', 'S', 'G', 'Q', 'T', 'N', 'H', 'E', 'D', 'R', 'K', 'C']


def read_domain_file(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Read a SHEPHARD-style IDR table (accession, start, end, type; 1-indexed inclusive).

    Parameters
    ----------
    path : pathlib.Path
        Table written by ``build_idr_tables.py``.

    Returns
    -------
    dict of str to list of tuple
        ``{accession: [(start, end), ...]}`` for IDR rows.
    """
    if not path.exists():
        raise FileNotFoundError(f"MISSING: {path}. Run `python build_idr_tables.py` first.")
    idrs = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            accession, start, end, domain_type = line.rstrip('\n').split('\t')
            assert domain_type == 'IDR', f'unexpected domain type {domain_type} in {path}'
            idrs[accession].append((int(start), int(end)))
    return dict(idrs)


def load_set(name: str) -> tuple[dict[str, str], dict[str, list[tuple[int, int]]]]:
    """Load sequences and IDRs for ``'proteome'`` or ``'tfs'``.

    Returns
    -------
    tuple
        ``(sequences, idrs)``; every IDR accession is checked to be in ``sequences``.
    """
    if name == 'proteome':
        # a few human proteome entries contain non-standard residues (e.g. X), converted
        sequences = read_uniprot_fasta(HUMAN_PROTEOME_FASTA, invalid_sequence_action='convert')
        idrs = read_domain_file(HUMAN_PROTEOME_IDRS)
    elif name == 'tfs':
        sequences = read_uniprot_fasta(HUMAN_TF_FASTA, invalid_sequence_action='fail')
        idrs = read_domain_file(HUMAN_TF_IDRS)
    else:
        raise ValueError(f"name must be 'proteome' or 'tfs', not {name!r}")
    missing = set(idrs) - set(sequences)
    assert not missing, f'{len(missing)} IDR accessions not in the FASTA, e.g. {sorted(missing)[:3]}'
    return sequences, idrs


def idr_residue_counts(sequences: dict[str, str], idrs: dict[str, list[tuple[int, int]]]) -> dict[str, int]:
    """Number of residues inside IDRs for every protein (0 if it has none)."""
    return {acc: sum(end - start + 1 for start, end in idrs.get(acc, [])) for acc in sequences}


def composition(sequences: dict[str, str], idrs: dict[str, list[tuple[int, int]]] | None = None) -> dict[str, float]:
    """Amino-acid composition of whole proteins, or of their IDRs only.

    Parameters
    ----------
    sequences : dict of str to str
        ``{accession: sequence}``.
    idrs : dict or None, optional
        If given, only residues inside these IDRs are counted.

    Returns
    -------
    dict of str to float
        Fraction of each of the 20 amino acids (sums to 1).
    """
    counts = dict.fromkeys(AMINO_ACIDS, 0)
    for acc, seq in sequences.items():
        regions = [seq] if idrs is None else [seq[start - 1:end] for start, end in idrs.get(acc, [])]
        for region in regions:
            for aa in region:
                counts[aa] += 1
    total = np.sum(list(counts.values()))
    return {aa: counts[aa] / total for aa in AMINO_ACIDS}
