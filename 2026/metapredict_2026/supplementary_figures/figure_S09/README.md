# Figure S9 — a stringent PDB-disorder-monomer benchmark

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S09A.ipynb` | Top 40 of 66 CAID3 disorder predictors by AUC on the PDB-disorder-monomer set. |
| B | `figure_S09B.ipynb` | Top 40 by APS. |
| C | `figure_S09C.ipynb` | Top 40 by Fmax. |
| D | `figure_S09D.ipynb` | Top 40 by MCC at each predictor's Fmax cutoff. |
| A–D | `figure_S09.ipynb` | The four panels stacked with panel letters (`figures/figure_S09.pdf`). |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code, `caid3_panels.py` loads the data and draws the panels, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (~15 s each). Each panel notebook writes `figures/figure_S09<panel>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid3/`, `metapredict_caid/` and `pdb_cache/monomer_labels.json` from `../../data/` (see its README).

**The benchmark.** A residue is ordered if it is modelled in at least one PDB entry of the isolated, monomeric protein, and disordered if it lies inside a monomeric construct but is modelled in none. All other residues are excluded, leaving 134 proteins with 32,904 ordered and 3,807 disordered residues (8.6 : 1). The labels come from the PDB alone (built by `../../data/pdb_cache_pipeline/run_pdb_coverage.py`), so 835 of the disordered residues are modelled in a homo-oligomer or complex. Scoring is as for Figure S3.

**Result** (rank of 66 for AUC, APS, Fmax, MCC): V3 0.791 (11), 0.299 (12), 0.412 (9), 0.340 (10); V2 0.759 (19), 0.255 (21), 0.359 (18), 0.275 (20); V1 0.693 (43), 0.206 (42), 0.284 (46), 0.185 (46). AlphaFold-pLDDT is first on all four. One protein, DP04052 (3,056 residues), holds 1,257 of the 3,807 disordered residues, and 7 of the 66 predictors did not predict it.
