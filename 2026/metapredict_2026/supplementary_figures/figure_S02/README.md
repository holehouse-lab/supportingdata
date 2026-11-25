# Figure S2 — building the combined scores used to train metapredict V3

A worked example on Chromogranin-A (P05059): AF2 pLDDT is rescaled and inverted, combined with the metapredict V1 score to give the combined score that V3 was trained on, and shown next to V3's own prediction.

| File | Contents |
| --- | --- |
| `figure_S02.ipynb` | Writes the three plot components: `figures/figure_S02_inputs.pdf` (pLDDT, rescaled pLDDT, V1), `figures/figure_S02_combined.pdf` and `figures/figure_S02_v3.pdf`. |
| `hybrid_scores.py` | `meta_predict_hybrid_v3` and `moving_average`, copied unchanged from `metapredict/analysis/creating_V2_and_V3_scores.py` in the metapredict repository. |
| `figtools.py` | Shared style. |

**Running:** open the notebook from this directory and run all cells. Requires metapredict 3.1.x, numpy, scipy, tqdm and matplotlib. V1 and V3 scores are predicted live; pLDDT is read from `../../data/af2_plddt/P05059_af2_plddt.tsv` (see the README there).

**How the combined score is built (Methods)**

1. pLDDT/100 is rescaled with a base of 0.35 and a top of 0.95, (pLDDT/100 − 0.35)/(0.95 − 0.35), then inverted so high values mean disorder.
2. At each residue, if either the rescaled pLDDT or the V1 score is above 0.5 the higher is kept; otherwise the lower is kept.
3. A Savitzky–Golay filter (window 7, order 3) is applied, values are clipped to 0–1, and a 25-residue moving average is applied.

This was applied to 455,666 Swiss-Prot sequences with AF2 models to make the V3 training set.

**Notebook output vs the assembled figure.** They differ. The SI figure is the original 2024 figure: its plots have blue pLDDT traces, guide lines at 0.2–0.8 and plot titles, and its BRNN schematic was drawn by hand. The notebook draws the same three plots in the manuscript style (pLDDT in black, V1 orange, V3 crimson, only the 0.5 line), with V3 predicted by metapredict 3.1.0, and does not draw the schematic.
