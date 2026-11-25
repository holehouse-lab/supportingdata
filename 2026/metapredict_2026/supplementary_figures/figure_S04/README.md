# Figure S4 — metapredict V1, V2 and V3 against the top CAID3 predictors on five proteins

`figure_S04.ipynb` draws per-residue disorder scores for DP04023, DP04296, DP04242, DP04253 and DP04340 from metapredict V1, V2 and V3 and the five CAID3 predictors with the highest AUC (PUNCH2, PUNCH2-Light, AlphaFold-rsa, SPOT-Disorder2, AlphaFold3-rsa). It writes `figures/figure_S04.pdf`.

The bar under each trace gives the Disorder-PDB labels: disordered (orange), ordered (green) and unscored (hatched grey; residues 1–9 of DP04023).

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code, `caid3_panels.py` loads the data and draws the traces (the same `trace_panel` as Figure 2C), and `figtools.py` is the shared style.

**Running:** open the notebook from this directory and run all cells (~15 s). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid_references/`, `caid3/` and `metapredict_caid/` from `../../data/` (see its README).

**The five proteins** have 150–300 residues, are 21–67% disordered and have a metapredict V3 per-protein AUC of at least 0.97. They were chosen to show a range of disorder content and IDR architectures.
