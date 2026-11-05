# Figure S16 — disorder in human transcription factors

**A:** each human transcription factor (Lambert et al. 2018) drawn as its IDR (red) and non-IDR (black) residues, ordered by length. The inset compares the fraction disordered per protein for the human proteome and for transcription factors. **B:** amino-acid composition of the proteome, proteome IDRs, transcription factors and transcription-factor IDRs.

| File | Contents |
| --- | --- |
| `build_idr_tables.py` | **Run first.** Predicts IDRs with metapredict 3.1.x (V3, default settings, CPU) for the proteome and TF FASTA files and writes `HP_IDRs.tsv` and `TF_IDRs.tsv`, all in `../../data/human_transcription_factors/`. Takes about 1 min. |
| `figure_S16A.ipynb` | Panel A: `figures/figure_S16A.pdf` and the inset `figures/figure_S16A_inset.pdf`. |
| `figure_S16B.ipynb` | Panel B: `figures/figure_S16B.pdf`. |
| `tf_idrs.py` | Loads sequences and IDR tables; per-protein disorder and composition. |
| `optional_shephard_annotation.py` | Optional: the original 2024 SHEPHARD route to the same tables. Not needed for the figure. |
| `figtools.py` | Shared manuscript style (copy of `shared/figtools.py`). |

**Running:** `python build_idr_tables.py`, then run the notebooks from this directory (metapredict 3.1.x, protfasta, numpy, matplotlib). The FASTA inputs and IDR tables are in `../../data/human_transcription_factors/` (see its README). Fraction disordered for a protein is the number of residues inside its IDRs divided by its length.

**Changes vs the published panel (2024 predictions)**
- **Numbers are unchanged.** The IDR tables from 3.1.0 match the 2024 ones except for one boundary (Q13277: IDR 1–20 is now 1–19). The medians are unchanged: 0.236 for the proteome and 0.569 for TFs.
- **The inset marks the median.** It is labelled with the median but used to draw the mean line; it now draws the median. The proteome is black and TFs are blue, matching B.
- **B uses the xmgrace colours.** The four sets are black, red, blue and green, all solid.
- **No SHEPHARD.** The notebooks use only metapredict output. The TF set is the 1,627 reviewed UniProt entries mapped from Lambert et al.'s gene names.
