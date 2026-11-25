# Figure 3C — viruses that infect prokaryotes vs eukaryotes

Each point is one viral proteome, grouped by host: prokaryote-infecting (archaeal + bacterial hosts), *E. coli*-infecting, eukaryote-infecting and *H. sapiens*-infecting; black bars mark the mean. This was Figure 2C in the original analysis.

| File | Contents |
| --- | --- |
| `Fig2C.ipynb` | Draws the panel and runs the statistics quoted in the Methods. |
| `host_to_tax_id.py` | Splits viruses by host using Virus-Host DB and writes the five per-host tables. |
| `prokaryote_vs_eukaryote_infecting_viruses_AFPP.pdf` | The panel. |
| `fig_infecting_different_orgs.ipynb` | An earlier, wider exploration (extra host groups). It needs five per-host tables that were never distributed, so it stops with `FileNotFoundError`; kept as a record only. |

**Running**

1. `python host_to_tax_id.py` (optional; its outputs are provided). It reads `virushostdb.tsv`, `org_id_to_proteome_id.tsv` and the per-proteome disorder table, buckets each virus by the first rank of its host lineage (Eukaryota, Archaea, Bacteria), flags hosts with taxonomy ID 9606 (human) and 562 (*E. coli*), and writes `{archaeal,bacterial,eukaryotic,human,ecoli}_infecting_viruses.tsv`.
2. Open `Fig2C.ipynb` from this directory and run all cells (numpy, scipy, scikit-learn, matplotlib, seaborn and metapredict are imported; no prediction is run). It is set to `disorder_vals = 'per_protein'` (AFPP) and writes the PDF here.

All tables are in `../../../data/viral_proteomes/` (`DATA_DIR` in the notebooks; see the README there).

**Statistics (Methods, "Statistical testing"), on AFPP**

| Comparison | n | KS D (P) | Mann–Whitney P | CLES |
| --- | --- | --- | --- | --- |
| Prokaryote- vs eukaryote-infecting | 4,108 vs 2,666 | 0.311 (5 × 10⁻¹³⁹) | 2 × 10⁻⁸⁶ | 0.6415 |
| *E. coli*- vs *H. sapiens*-infecting | 206 vs 211 | 0.582 (4 × 10⁻³³) | 4 × 10⁻²² | 0.7742 |

Shapiro–Wilk rejects normality for all four groups (P < 10⁻⁴). The human and *E. coli* groups are subsets of the eukaryote and prokaryote groups.

**Notebook output vs the assembled figure.** The notebook draws the host groups along the x axis. In the assembled figure the panel is rotated and its axis title and group labels retyped.

**Fixed in the co-author code**

- The plotting cell was left on `'total_proteome'` (TFAP) although the figure and statistics use AFPP; it is now `'per_protein'`.
- The jitter was unseeded; `np.random.seed(42)` now precedes the plotting loop in both notebooks (cosmetic only).
- `host_to_tax_id.py` checked for 9 fields per row but read the tenth; it now finds the columns from the header and checks the row length against them.
