# Figure S13 — Figure 3 recomputed with total fraction across proteome (TFAP)

This is the TFAP version of Figure 3: disordered residues divided by all residues in each proteome, instead of the average per-protein fraction (AFPP).

> **Data predate metapredict 3.1.0.** Panels B and C read per-proteome disorder from the 2024 metapredict V3 prediction of 23,129 UniProt reference proteomes. That run is not repeated here, so these values come from the earlier metapredict release.

| Path | Contents |
| --- | --- |
| `panel_A_tree/` | Panel A (family-level tree): the published panel and the R script, copied from the original analysis; its per-family table and tree are in `../../data/tree_of_life/`. |
| `figure_S13B.ipynb` | Panel B: viral lineages (≥10 proteomes; single-gene proteomes excluded). |
| `figure_S13C.ipynb` | Panel C: viruses grouped by host (prokaryote, *E. coli*, eukaryote, *H. sapiens*). |
| `viral_disorder.py` | Table readers and jitter/strip-plot functions taken from the original notebook. |
| `figtools.py` | Shared manuscript style (copy of `shared/figtools.py`). |
| `host_to_tax_id.py` | Original script that splits viruses by host domain (Eukaryota / Archaea / Bacteria) and picks out *H. sapiens* (taxid 9606) and *E. coli* (taxid 562) hosts; writes the five per-host tables to `../../data/viral_proteomes/`. |

**Running:** open a notebook from this directory and run all cells (numpy, scipy, matplotlib; no prediction is run). Jitter is seeded (`JITTER_SEED`). The input tables are in `../../data/viral_proteomes/` (shared with Figure 3B–C; see the README there).

**Differences from the published panel**
- **B now actually plots TFAP.** The original notebook set the metric to TFAP but never passed it on, so the grouping fell back to AFPP and the published S6B duplicated Figure 3B (Figure 2B in the published version). TFAP means are lower; for example, Anelloviridae goes from 0.626 to 0.510. The group order is unchanged.
- **B drops the non-taxonomic bins in code.** The "unclassified …" and "environmental samples" bins were removed by hand in the published panel; here the code removes them.
- **C has tick labels.** Group names are drawn as x tick labels instead of being added at assembly.
- **Colours are unchanged.** The default matplotlib colour cycle is kept for the 9 lineages and 4 host groups (more than the 5 xmgrace colours), so the colours match Figure 3.
- **Orientation is unchanged.** B and C are drawn vertically and are rotated 90° when the figure is assembled, as before.
