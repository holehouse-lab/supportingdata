"""Regenerate the human proteome and human transcription factor IDR tables.

Predicts IDRs with metapredict (V3 network, default domain-decomposition settings)
for every protein in the human reference proteome and in the human transcription
factor set, and writes them as SHEPHARD-style domain files (one IDR per line:
UniProt accession, start, end, ``IDR``; 1-indexed, inclusive). These are the
inputs read by the Figure S16 notebooks.

Run from this directory (takes ~1-2 min on a laptop CPU)::

    python build_idr_tables.py

The defaults match the IDR annotation originally done with SHEPHARD's
``metapredict_api.annotate_proteome_with_disordered_domains`` (V3, the network's
own disorder threshold, minimum IDR 12, minimum folded domain 50, gap closure 10),
so the output format and content are interchangeable with the 2024 tables.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import metapredict as meta
import protfasta

from figtools import SHARED_DATA_DIR, require_metapredict_3_1

# Inputs and the IDR tables written here live in data/human_transcription_factors/
# at the top of the repository (see the README there for provenance).
DATA_DIR = SHARED_DATA_DIR / "human_transcription_factors"
HUMAN_PROTEOME_FASTA = DATA_DIR / "UP000005640_9606.fasta"
HUMAN_TF_FASTA = DATA_DIR / "human_TFs_full_2024_10_29.fasta"
HUMAN_PROTEOME_IDRS = DATA_DIR / "HP_IDRs.tsv"
HUMAN_TF_IDRS = DATA_DIR / "TF_IDRs.tsv"

# CPU by default so the tables are bit-for-bit reproducible across machines
# (GPU/MPS predictions differ from CPU at the ~1e-6 level).
DEFAULT_DEVICE = "cpu"


def read_uniprot_fasta(path: Path, invalid_sequence_action: str) -> dict[str, str]:
    """Read a UniProt FASTA file keyed by accession.

    Parameters
    ----------
    path : pathlib.Path
        UniProt FASTA file with headers like ``>sp|P33748|MSN2_YEAST ...``.
    invalid_sequence_action : str
        Passed to ``protfasta.read_fasta``; ``'convert'`` maps non-standard
        residues (e.g. ``X``) onto standard ones, as in the original analysis.

    Returns
    -------
    dict of str to str
        ``{accession: sequence}``.
    """
    records = protfasta.read_fasta(str(path), invalid_sequence_action=invalid_sequence_action)
    sequences = {header.split("|")[1]: seq for header, seq in records.items()}
    assert len(sequences) == len(records), f"duplicate accessions in {path}"
    return sequences


def predict_idrs(sequences: dict[str, str], device: str) -> dict[str, list[tuple[int, int]]]:
    """Predict IDR boundaries for every sequence.

    Parameters
    ----------
    sequences : dict of str to str
        ``{accession: sequence}``.
    device : str
        Torch device for prediction (``'cpu'``, ``'mps'`` or ``'cuda'``).

    Returns
    -------
    dict of str to list of tuple
        ``{accession: [(start, end), ...]}`` with 1-indexed, inclusive positions.
    """
    predictions = meta.predict_disorder(sequences, version="V3", return_domains=True,
                                        device=device, show_progress_bar=True)
    idrs = {}
    for accession, disorder_object in predictions.items():
        # metapredict boundaries are 0-indexed, half-open [start, end)
        idrs[accession] = [(start + 1, end) for start, end in disorder_object.disordered_domain_boundaries]
    return idrs


def write_domain_file(idrs: dict[str, list[tuple[int, int]]], path: Path) -> None:
    """Write IDRs as a SHEPHARD domain file (accession, start, end, IDR).

    Parameters
    ----------
    idrs : dict of str to list of tuple
        Output of :func:`predict_idrs`.
    path : pathlib.Path
        Output TSV path.
    """
    n_idrs = 0
    with open(path, "w") as fh:
        for accession, boundaries in idrs.items():
            for start, end in boundaries:
                fh.write(f"{accession}\t{start}\t{end}\tIDR\n")
                n_idrs += 1
    print(f"Wrote {path}: {len(idrs)} proteins, {n_idrs} IDRs")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--device", default=DEFAULT_DEVICE, help="cpu, mps or cuda")
    args = parser.parse_args()
    require_metapredict_3_1()

    # a handful of human proteome entries contain non-standard residues (e.g. X)
    proteome = read_uniprot_fasta(HUMAN_PROTEOME_FASTA, invalid_sequence_action="convert")
    write_domain_file(predict_idrs(proteome, args.device), HUMAN_PROTEOME_IDRS)

    tfs = read_uniprot_fasta(HUMAN_TF_FASTA, invalid_sequence_action="fail")
    write_domain_file(predict_idrs(tfs, args.device), HUMAN_TF_IDRS)


if __name__ == "__main__":
    main()
