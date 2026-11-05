# Figure S2 — building the combined scores used to train metapredict V3

A worked example on Chromogranin-A (P05059). AF2 pLDDT is rescaled and inverted, combined with metapredict V1 disorder to give the "combined score" that V3 was trained on, and the result is shown next to V3's own prediction. The BRNN schematic in the middle of the published figure was drawn by hand and is not generated here.

| File | Contents |
| --- | --- |
| `figure_S02.ipynb` | Writes the three plot components: `figures/figure_S02_inputs.pdf` (pLDDT, rescaled pLDDT, V1), `figures/figure_S02_combined.pdf` and `figures/figure_S02_v3.pdf`. |
| `hybrid_scores.py` | `meta_predict_hybrid_v3` and `moving_average`, copied verbatim from `metapredict/analysis/creating_V2_and_V3_scores.py` (no longer shipped in the installed package). |
| `figtools.py` | Shared manuscript style/colours (copy of `shared/figtools.py`). |

**Running:** open the notebook from this directory and run all cells. Requires metapredict 3.1.x, numpy, scipy, tqdm and matplotlib. V1 and V3 scores are predicted live; the pLDDT values are read from `../../data/af2_plddt/P05059_af2_plddt.tsv` (see the README there).

**How the combined score is built (Methods)**

1. pLDDT / 100 is rescaled with base 0.35 and top 0.95, so (pLDDT/100 − 0.35)/(0.95 − 0.35), and then inverted (1 − value), so that high values mean disorder.
2. At each residue, if either the rescaled pLDDT or the metapredict V1 score is above 0.5 the larger one is kept; otherwise the smaller one is kept.
3. A Savitzky–Golay filter (window 7, order 3) is applied, values are clipped to [0, 1], and a 25-residue moving average is applied.

These steps are `meta_predict_hybrid_v3` in `hybrid_scores.py`. It was applied to 455,666 de-duplicated Swiss-Prot sequences with AF2 models to make the V3 training set.

**Differences from the published panel**
- Only the 0.5 line is drawn (dashed); the 0.2/0.4/0.6/0.8 guide lines are gone, and so is the 50 line on the raw pLDDT plot.
- pLDDT traces are black instead of blue (blue is now the V2 colour); V1 is orange and V3 crimson as everywhere else.
- Plot titles ("Combined Score", "Disorder Score") are now y-axis labels.
