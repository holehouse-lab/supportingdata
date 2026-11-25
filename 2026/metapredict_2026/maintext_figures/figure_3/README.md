# Figure 3 — disorder across the tree of life

This figure was made by co-authors. The code is copied from the original analysis, where it was Figure 2, so the files inside keep their original names (`Figure_2B/`, `Fig2C.ipynb`, ...). The only changes are the data paths and the bug fixes listed in each subdirectory's README.

| Panel | Directory | Contents |
| --- | --- | --- |
| A | `Figure_2A_just_analysis_no_notebook/` | R analysis (taxize) that makes the per-family disorder means and the family tree, which was drawn in iTOL. No plotting notebook. |
| B | `Figure_2B/` | Viral proteome disorder by viral lineage. |
| C | `Figure_2C/` | Viruses that infect prokaryotes vs eukaryotes, with the statistics. |
| — | `tree_assets/` | The panel A tree without lettering (PDF, PNG, SVG). |

**Data.** Inputs are in `../../data/tree_of_life/` (A) and `../../data/viral_proteomes/` (B, C), shared with Figure S13. The disorder values are the 2024 metapredict V3 predictions of 23,129 UniProt reference proteomes (10,561 viral, 12,568 cellular); they are not regenerated here. Figure 3 uses the average fraction disordered per protein (AFPP); Figure S13 uses the total fraction across the proteome (TFAP).

The notebooks in `Figure_2B/` and `Figure_2C/` keep their original 2024 style and colours, and write their PDFs next to the notebook instead of into `figures/`.

**Notebook output vs the assembled figure.** No code draws panel A: the tree was drawn in iTOL and the organism images added by hand. The B and C notebooks draw vertical strip plots with the groups along the x axis; in the assembled figure they are rotated, and the axis title ("Fraction of proteome disordered") and group labels were retyped in Illustrator. The values and means are the same.
