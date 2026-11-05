# Figure S11 — CAID3 Disorder-PDB scored by IDR length

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S11A.ipynb` | AUC of all 66 CAID3 disorder predictors in four IDR-length bins; the top five CAID3 methods and metapredict V1/V2/V3 as lines, metapredict V3's rank above each bin. The legend is drawn here only. |
| B | `figure_S11B.ipynb` | APS, as in A. |
| C | `figure_S11C.ipynb` | Fmax, as in A. |
| D | `figure_S11D.ipynb` | MCC at each method's Fmax threshold, as in A. |
| E | `figure_S11E.ipynb` | Top 40 by MCC, IDRs of 50–99 residues; an asterisk marks a predictor significantly ahead of metapredict V3 (paired bootstrap). |
| F | `figure_S11F.ipynb` | As E, IDRs of at least 100 residues. |
| A–F | `figure_S11.ipynb` | All six panels in one figure with panel letters (`figures/figure_S11.pdf`): A–D in a 2 × 2 grid with the legend beside B, E and F below at full width. |

`idr_length.py` makes the bins, runs the paired bootstrap and draws the panels. `caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (A–D ~1 min each, E–F ~20 s). Reads the CAID3 reference set (`caid_references/`, from which `idr_length.py` makes the IDR-length bins), the CAID3 predictions (`caid3/`) and the metapredict 3.1.0 predictions (`metapredict_caid/`) from the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Each panel notebook writes `figures/figure_S11<panel>.pdf`; `figure_S11.ipynb` (~2 min) draws all six into `figures/figure_S11.pdf`.

**The bins:** an IDR is a maximal run of consecutive residues labelled disordered in the reference. Each bin keeps all 67,838 ordered residues and only the disordered residues in IDRs of that length; every other disordered residue is excluded. Each bin is then scored exactly as Figure S3 scores the full set (pooled per-residue AUC, APS, Fmax, and MCC at each method's Fmax threshold, chosen within the bin).

| IDR length (residues) | IDRs | Disordered residues |
| --- | --- | --- |
| <20 | 158 | 2,123 |
| 20–49 | 135 | 4,042 |
| 50–99 | 68 | 4,801 |
| ≥100 | 94 | 20,435 |

**Significance (E–F):** targets are resampled with replacement 1,000 times (seed 20261030); every predictor ahead of V3 and V3 itself are scored on the same resampled targets (those the predictor covers), and a predictor is significantly ahead if the central 95% interval of (predictor − V3) lies above zero. The bootstrap uses a weighted re-implementation of the pooled metrics (`idr_length.weighted_metrics`), which the notebooks check against the CAID3 scoring code before running.

**Result with metapredict 3.1.0** (V3 rank of 66, AUC / APS / Fmax / MCC): <20 aa 49 / 51 / 47 / 47; 20–49 aa 33 / 35 / 32 / 30; 50–99 aa 13 / 6 / 5 / 5; ≥100 aa 7 / 5 / 4 / 4. For IDRs of 50–99 residues no method is significantly ahead of V3 on APS, Fmax or MCC (AlphaFold-rsa and AlphaFold3-rsa are on AUC). For IDRs of at least 100 residues only PUNCH2 is significantly ahead on APS, Fmax and MCC (PUNCH2 and PUNCH2-Light on AUC); MCC PUNCH2 0.866, V3 0.835.

New figure; not part of the original analysis.
