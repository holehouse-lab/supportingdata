# Figure S13 — Figure 3 with the total fraction across the proteome

Figure 3 re-plotted with each proteome's disorder measured as disordered residues divided by all residues (total fraction across the proteome, TFAP), instead of the average fraction per protein (AFPP). The data are the 2024 metapredict V3 predictions of the UniProt reference proteomes, which are not regenerated here.

| Panel | File | Contents |
| --- | --- | --- |
| A | `panel_A_tree/` | The family-level tree (1,480 families) and the R script behind it, copied from the Figure 3A analysis. |
| B | `figure_S13B.ipynb` | Viral proteomes by lineage (lineages with at least 10 proteomes; single-protein proteomes excluded). |
| C | `figure_S13C.ipynb` | Viral proteomes by host, and the two KS tests quoted in the legend. |

`viral_disorder.py` holds the table readers and strip-plot functions from the original notebooks, `figtools.py` is the shared style, and `host_to_tax_id.py` is the original script that writes the five per-host tables (their outputs are provided).

**Running:** open a notebook from this directory and run all cells (numpy, scipy and matplotlib; no prediction is run). Each writes `figures/figure_S13<panel>.pdf`. The tables are in `../../data/viral_proteomes/` (shared with Figure 3B–C; see the README there). Jitter is seeded (`JITTER_SEED`).

**Notebook output vs the assembled figure.** No code draws panel A: the tree was drawn in iTOL and finished by hand. The B and C notebooks draw vertical strip plots with the groups along the x axis; in the assembled figure they are rotated. The points and means are the same.

**Result**

- B: TFAP means run from 0.047 (Adnaviria) to 0.510 (Anelloviridae); the lineage order is the same as with AFPP in Figure 3B (0.069 to 0.626). The notebook prints both.
- C: mean TFAP is 0.111 for prokaryote-infecting viruses (n = 4,108), 0.098 for *E. coli*-infecting (206), 0.183 for eukaryote-infecting (2,666) and 0.271 for *H. sapiens*-infecting (211). KS tests: prokaryote vs eukaryote D = 0.293, P = 1 × 10⁻¹²²; *E. coli* vs *H. sapiens* D = 0.610, P = 1 × 10⁻³⁶.

**Note on B.** In the original 2024 analysis the TFAP setting changed only the axis label, so that panel showed AFPP values and duplicated Figure 3B (see `maintext_figures/figure_3/Figure_2B/README.md`). `figure_S13B.ipynb` plots TFAP, and drops the non-taxonomic bins ("unclassified …", "environmental samples") in code.
