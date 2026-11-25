# Figure S11 — CAID3 Disorder-PDB scored by IDR length

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S11A.ipynb` | AUC of all 66 CAID3 disorder predictors in four IDR-length bins, with metapredict V3's rank above each bin. Carries the legend. |
| B | `figure_S11B.ipynb` | APS, as in A. |
| C | `figure_S11C.ipynb` | Fmax, as in A. |
| D | `figure_S11D.ipynb` | MCC at each predictor's Fmax cutoff, as in A. |
| E | `figure_S11E.ipynb` | Top 40 by MCC for IDRs of 50–99 residues; an asterisk marks a predictor significantly more accurate than metapredict V3. |
| F | `figure_S11F.ipynb` | As E, for IDRs of 100 residues or more. |
| A–F | `figure_S11.ipynb` | All six panels with panel letters (`figures/figure_S11.pdf`). |

`idr_length.py` makes the bins, runs the bootstrap and draws the panels. `caid3_panels.py`, `caid_eval.py` and `predictor_classes.py` are the CAID3 code, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (A–D ~1 min each, E–F ~20 s, the combined figure ~2 min). Each panel notebook writes `figures/figure_S11<panel>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid_references/`, `caid3/` and `metapredict_caid/` from `../../data/` (see its README).

**The bins (Methods).** An IDR is a run of residues labelled disordered in CAID3 Disorder-PDB; unscored residues end a run, and the shortest IDR is 10 residues. Each bin keeps all 67,838 ordered residues and only the disordered residues of IDRs in that bin, and is scored as in Figure S3 with the Fmax cutoff chosen within the bin.

| IDR length (residues) | IDRs | Disordered residues | V3 rank of 66 (AUC / APS / Fmax / MCC) |
| --- | --- | --- | --- |
| 10–19 (axis label `<20`) | 158 | 2,123 | 49 / 51 / 47 / 47 |
| 20–49 | 135 | 4,042 | 33 / 35 / 32 / 30 |
| 50–99 | 68 | 4,801 | 13 / 6 / 5 / 5 |
| ≥100 | 94 | 20,435 | 7 / 5 / 4 / 4 |

**Significance (E, F).** A paired bootstrap over proteins: 1,000 resamples of the 319 proteins with replacement (seed 20261030), the same resamples for every predictor, each predictor and V3 scored on the drawn proteins the predictor covers. A predictor is marked if the 2.5th percentile of its difference from V3 is above zero. By MCC, no predictor is significantly ahead of V3 for IDRs of 50–99 residues, and only PUNCH2 is for IDRs of 100 residues or more (MCC 0.866 vs 0.835). The bootstrap uses a fast weighted version of the pooled metrics (`idr_length.weighted_metrics`), which the notebooks check against the CAID3 scoring code first.
