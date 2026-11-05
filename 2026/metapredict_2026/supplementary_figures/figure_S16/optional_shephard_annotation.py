"""OPTIONAL: the original (2024) SHEPHARD route to the Figure S16 IDR tables.

Not needed to make the figure - ``build_idr_tables.py`` produces the same tables
with metapredict alone. Kept so the original analysis path is documented. Requires
``pip install shephard`` (https://shephard.readthedocs.io).

All files are in ``data/human_transcription_factors/`` at the top of the repository
(``DATA_DIR`` below). Input provenance (manual UniProt steps, not scripted):

* ``UP000005640_9606.fasta`` - human reference proteome, one sequence per gene,
  downloaded from https://www.uniprot.org/proteomes/UP000005640 on 2024-10-29.
* ``human_TFs_full_2024_10_29.fasta`` - reviewed human UniProt entries for the
  gene names in ``human_TFs_gene_name.txt`` (Lambert et al. 2018, Cell 172,
  650-665), retrieved with the UniProt ID Mapping tool
  (https://www.uniprot.org/id-mapping) on 2024-10-29.

Run from this directory: ``python optional_shephard_annotation.py``. Writes
``HP_IDRs_shephard.tsv`` and ``TF_IDRs_shephard.tsv`` (same format as
``build_idr_tables.py``) and ``human_tf_uniprot_IDs.txt``.
"""

from shephard import interfaces
from shephard.apis import metapredict_api, uniprot

from build_idr_tables import DATA_DIR
from figtools import require_metapredict_3_1


def main() -> None:
    """Annotate IDRs in the TF set and the human proteome with SHEPHARD and write them out."""
    require_metapredict_3_1()

    # Read the TFs and the human proteome as SHEPHARD Proteome objects; the
    # proteome has a few 'X' residues, which are converted.
    TFs = uniprot.uniprot_fasta_to_proteome(str(DATA_DIR / 'human_TFs_full_2024_10_29.fasta'))
    HP = uniprot.uniprot_fasta_to_proteome(str(DATA_DIR / 'UP000005640_9606.fasta'),
                                           invalid_sequence_action='convert')

    # annotate IDRs with metapredict (V3, default settings)
    metapredict_api.annotate_proteome_with_disordered_domains(TFs)
    metapredict_api.annotate_proteome_with_disordered_domains(HP)

    interfaces.si_domains.write_domains(TFs, str(DATA_DIR / 'TF_IDRs_shephard.tsv'))
    interfaces.si_domains.write_domains(HP, str(DATA_DIR / 'HP_IDRs_shephard.tsv'))

    # save human TF UniProt IDs
    with open(DATA_DIR / 'human_tf_uniprot_IDs.txt', 'w') as fh:
        for protein in TFs.proteins:
            fh.write(f"{protein}\n")


if __name__ == '__main__':
    main()
