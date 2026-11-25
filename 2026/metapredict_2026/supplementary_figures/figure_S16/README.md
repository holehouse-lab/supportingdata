# Figure S16 — disorder in human transcription factors

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S16A.ipynb` | Each transcription factor as its residues inside IDRs (red) and outside them (black), ordered by length (`figures/figure_S16A.pdf`), and the inset comparing the fraction disordered per protein in the human proteome and in transcription factors (`figures/figure_S16A_inset.pdf`). |
| B | `figure_S16B.ipynb` | Amino-acid composition of the proteome, proteome IDRs, transcription factors and transcription-factor IDRs (`figures/figure_S16B.pdf`). |

| File | Contents |
| --- | --- |
| `build_idr_tables.py` | Predicts IDRs with metapredict 3.1.x (V3, default settings, CPU) for the proteome and transcription-factor FASTA files and writes `HP_IDRs.tsv` and `TF_IDRs.tsv`. Takes about 1 min. |
| `tf_idrs.py` | Loads sequences and IDR tables; per-protein disorder and composition. |
| `optional_shephard_annotation.py` | The original 2024 SHEPHARD route to the same tables. Not needed for the figure. |
| `figtools.py` | Shared style. |

**Running:** open a notebook from this directory and run all cells (metapredict 3.1.x, protfasta, numpy, matplotlib). The FASTA files and IDR tables are in `../../data/human_transcription_factors/` (see its README); the tables are provided, and `python build_idr_tables.py` rebuilds them.

**The transcription-factor set.** The list from Lambert *et al.* (2018) has 1,608 gene names. Mapping them to reviewed human UniProt entries (29 October 2024) gives 1,627 sequences, and the figure plots all 1,627.

**Result:** the fraction disordered of a protein is its residues inside IDRs divided by its length. The median is 0.569 for transcription factors and 0.236 for the human proteome (20,656 proteins), a ratio of 2.4.
