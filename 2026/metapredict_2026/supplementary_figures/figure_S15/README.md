# Figure S15 — metapredict V3 on synthetic and designed proteins

metapredict V3 predictions for 63 non-natural sequences from four studies; 61 of 63 are called correctly.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A, B, C | `figure_S15A.ipynb`, `figure_S15B.ipynb`, `figure_S15C.ipynb` | Disorder profiles of the hallucinated folded proteins 0738, 0515 and 0217 (Anishchenko *et al.* 2021). |
| D | `figure_S15D.ipynb` | Wild-type hnRNPA1 low-complexity domain and five designed variants (Pesce *et al.* 2024). |
| E | `figure_S15E.ipynb` | Fraction predicted disordered for 7 folded and 15 disordered random-sequence proteins (Tretyachenko *et al.* 2017). |
| F | `figure_S15F.ipynb` | Fraction predicted disordered for 32 designed disordered proteins (Hunter *et al.* 2026; SI reference 4). |

`synthetic_disorder.py` holds the plotting and fraction-disordered functions, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. Each writes `figures/figure_S15<panel>.pdf`. Requires metapredict 3.1.x, protfasta, numpy and matplotlib; all predictions are made live. The sequences for D–F are in `../../data/synthetic_proteins/` (see its README). The sequences for A–C are written in their notebooks as the expressed constructs; 0738 and 0217 start with a 20-residue His-tag.

**Scoring.** The fraction disordered is the fraction of residues inside a metapredict-defined disordered domain. A sequence counts as disordered if more than half of it is.

**Result:** A–C are predicted folded apart from the His-tags in A and C; all six sequences in D and all 32 in F are predicted disordered; in E all 7 folded and 13 of 15 disordered sequences are called correctly (4954 and 9121 are predicted mostly folded; MCC 0.821).

The dashed line in A and C marks residue 21, the first residue after the His-tag. The assembled figure uses the six panels as drawn.
