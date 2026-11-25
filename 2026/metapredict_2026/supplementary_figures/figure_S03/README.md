# Figure S3 — CAID3 Disorder-PDB rankings by AUC, APS, Fmax and MCC

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S03A.ipynb` | Top 40 of 66 CAID3 disorder predictors by AUC. |
| B | `figure_S03B.ipynb` | Top 40 by APS (average precision). |
| C | `figure_S03C.ipynb` | Top 40 by Fmax. |
| D | `figure_S03D.ipynb` | Top 40 by MCC at each predictor's Fmax cutoff. |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code, `caid3_panels.py` loads the data and draws the panels, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (~15 s each). Each writes `figures/figure_S03<panel>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid_references/`, `caid3/` and `metapredict_caid/` from `../../data/` (see its README).

**Scoring.** As for Figure 2A–B: the 63 CAID3 disorder predictors (binding and linker predictors excluded) plus metapredict V1–V3, pooled over the 31,401 disordered and 67,838 ordered residues of the 319 proteins. AUC and APS come from scikit-learn, Fmax is the highest F1 over all cutoffs, and MCC is taken at the Fmax cutoff, so it differs from the MCC of Figure 1C (each predictor's own binary calls).

**Result** (rank of 66 for AUC, APS, Fmax, MCC): V3 0.929 (14), 0.902 (7), 0.828 (8), 0.761 (7); V2 0.915 (30), 0.875 (18), 0.804 (19), 0.721 (18); V1 0.887 (45), 0.837 (38), 0.749 (43), 0.640 (41).
