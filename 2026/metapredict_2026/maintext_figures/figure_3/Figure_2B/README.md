# Figure 3B — disorder across viral lineages

Each point is one viral proteome, grouped by viral lineage (the first rank below "Viruses" in the UniProt lineage) and ordered by mean disorder; black bars mark the mean. This was Figure 2B in the original analysis.

| File | Contents |
| --- | --- |
| `Figure2B.ipynb` | Loads the per-proteome disorder values, groups them by lineage, filters and draws the strip plot. |
| `Fig2B_AFPP_min_10_per_group_no_single_gene_proteomes.pdf` | The panel. |

**Running:** open `Figure2B.ipynb` from this directory and run all cells (numpy, matplotlib, scikit-learn and metapredict are imported; no prediction is run). It reads three tables from `../../../data/viral_proteomes/` (`DATA_DIR` in the first code cell; see the README there) and writes the PDF here.

**Settings used for the figure**

- `disorder_vals = 'per_protein'`: the average fraction disordered per protein (AFPP).
- `remove_single_gene_proteome = True`: the 547 single-protein proteomes are dropped, leaving 10,014 of 10,561.
- `min_per_group = 10`: lineages with fewer than 10 proteomes are dropped.
- The non-taxonomic bins ("unclassified …" and "environmental samples") are left out, leaving nine lineages.

**Notes**

- `disorder_vals = 'total_proteome'` changes only the axis label and file name: the final cell calls `break_down_into_groups` without `vals`, so AFPP is always plotted. The TFAP version is Figure S13B (`supplementary_figures/figure_S13/`).
- Fixed here: the final cell used `pin_after` (the list of non-taxonomic bins) without defining it, so the notebook stopped with a `NameError` on a fresh kernel. The list is now defined in that cell.
- The point jitter is random and unseeded, so point positions change between runs; the values and means do not.
- The notebook draws the lineages along the x axis. In the assembled figure the panel is rotated and its axis title retyped.
