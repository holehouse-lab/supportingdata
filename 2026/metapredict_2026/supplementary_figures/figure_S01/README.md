# Figure S1 — metapredict V1 disorder vs AlphaFold2 pLDDT

metapredict V1 scores and AF2 pLDDT disagree about which regions are disordered, which motivated combining them. Panel A is EIF4EBP1 (Q13541) and panel B is Myozenin-1 (Q9NP98).

| File | Contents |
| --- | --- |
| `figure_S01A.ipynb`, `figure_S01B.ipynb` | One notebook per panel; each writes `figures/figure_S01<panel>.pdf`. |
| `disorder_vs_plddt.py` | The plotting function and the two thresholds: V1 scores above 0.42 and pLDDT below 70 are shaded. |
| `figtools.py` | Shared style. |

**Running:** open a notebook from this directory and run all cells. Requires metapredict 3.1.x, numpy and matplotlib. V1 scores are predicted live; pLDDT is read from `../../data/af2_plddt/` (see its README).
