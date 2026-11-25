# Figure S12 — the top CAID3 predictors differ by only a few residues per 100

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S12A.ipynb` | Correctly scored residues per 100 for the top 40 of 66 CAID3 disorder predictors, with a bracket from the best predictor to metapredict V3. |
| B | `figure_S12B.ipynb` | Extra residues scored incorrectly per 100 relative to the best predictor, with a bracket to the tenth-best. |

`caid_accuracy.py` reads the CAID files and scores them, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. Each writes `figures/figure_S12<panel>.pdf`. Requires metapredict 3.1.x, numpy and matplotlib. Reads `caid_references/` (reference set and list of disorder predictors), `caid3/` and `metapredict_caid/` from `../../data/` (see its README).

**Scoring.** Binary calls are as in Figure 1C: each predictor's submitted calls where available, otherwise a cutoff of 0.5 (0.581 and 0.312 for the AlphaFold2/AlphaFold3 RSA and pLDDT tracks); metapredict residues inside a DisorderedDomain are disordered. Accuracy is correct calls per 100 assessed residues, pooled over all 319 proteins.

**Result:** PUNCH2 is best at 91.5 correct per 100. metapredict V3 scores 89.6 (7th of 66, 1.8 fewer), V2 89.0 (9th) and V1 86.8 (23rd). The tenth-best predictor is 2.7 residues per 100 behind the best.
