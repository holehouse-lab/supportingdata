# Figure S3 — CAID3 Disorder-PDB rankings by AUC, APS, Fmax and MCC

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S03A.ipynb` | Top 40 CAID3 disorder predictors by AUC. |
| B | `figure_S03B.ipynb` | Top 40 by APS (average precision). |
| C | `figure_S03C.ipynb` | Top 40 by Fmax. |
| D | `figure_S03D.ipynb` | Top 40 by MCC at each method's Fmax threshold. |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data and draws the panels, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (~15 s each). Reads the CAID3 reference set (`caid_references/`), the CAID3 predictions (`caid3/`) and the metapredict 3.1.0 predictions (`metapredict_caid/`) from the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Each notebook writes `figures/figure_S03<panel>.pdf`.

**How the panels are scored.** Metrics are pooled over every covered, non-excluded residue of every target, reproducing the CAID3 Disorder-PDB table: AUC and APS from scikit-learn, Fmax as the maximum F1 over thresholds, MCC at the Fmax threshold. Binding and linker predictor flavours (`predictor_classes.py`) are excluded. The original re-analysis reproduced the published CAID3 top 10 by name and order, with AUC within 0.001.

**Result with metapredict 3.1.0** (rank of 66): V3 AUC 0.929 (14), APS 0.902 (7), Fmax 0.828 (8), MCC 0.761 (7); V2 0.915 (30), 0.875 (18), 0.804 (19), 0.721 (18); V1 0.887 (45), 0.837 (38), 0.749 (43), 0.640 (41).

**Differences from the published panel**
- metapredict is predicted with metapredict 3.1.0 instead of taken from the CAID3-distributed files (raw scores agree to ±0.001, ranks unchanged).
- One panel per metric, solid black bars, and no suptitle; the dataset name is written in the axes.
