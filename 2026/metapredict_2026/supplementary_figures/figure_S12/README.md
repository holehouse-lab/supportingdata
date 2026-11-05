# Figure S12 — the top CAID3 predictors differ by only a few residues per 100

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S12A.ipynb` | Correctly scored residues per 100 for the top 40 CAID3 disorder predictors; bracket from the best predictor to metapredict V3. |
| B | `figure_S12B.ipynb` | Extra residues scored incorrectly per 100 relative to the best predictor; bracket to the tenth-best. |

`caid_accuracy.py` holds the reader functions from the original S7 notebook and the CAID3 scoring, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. Reads the CAID3 reference set and predictor list (`caid_references/`), the CAID3 predictions (`caid3/`) and the metapredict 3.1.0 predictions (`metapredict_caid/`) from the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy and matplotlib. Each notebook writes `figures/figure_S12<panel>.pdf`.

**Scoring.** Each residue is called disordered or ordered from the method's own submitted binary calls; methods that submitted only raw scores are thresholded at 0.5 (AlphaFold RSA 0.581, AlphaFold pLDDT 0.312, the CAID2 thresholds), as in Figure 1C. Accuracy is (correct calls / assessed residues) × 100, pooled over all targets.

**Result with metapredict 3.1.0:** PUNCH2 91.5 correct per 100; metapredict V3 89.6 (7th of 66, 1.8 fewer than the best), V2 89.0 (9th), V1 86.8 (23rd); the tenth-best is 2.7 residues per 100 behind the best.

**Differences from the published panel**
- Computed on CAID3 Disorder-PDB instead of CAID2, top 40 shown, and V1/V2 are highlighted as well as V3.
- metapredict uses metapredict 3.1.0 predictions with 3.1.0's default binary calls (residues inside predicted IDRs are disordered).
- The original panel B bracketed the eleventh-ranked predictor while labelling it the tenth-best; this one uses the tenth-best.
