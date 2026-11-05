# Figure S9 — a stringent PDB-disorder-monomer benchmark

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S09A.ipynb` | Top 40 CAID3 disorder predictors by AUC on the PDB-disorder-monomer set. |
| B | `figure_S09B.ipynb` | Top 40 by APS. |
| C | `figure_S09C.ipynb` | Top 40 by Fmax. |
| D | `figure_S09D.ipynb` | Top 40 by MCC at each method's Fmax threshold. |
| A–D | `figure_S09.ipynb` | The four panels stacked into one figure with panel letters (`figures/figure_S09.pdf`). |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data and draws the panels, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (~15 s each). Reads the CAID3 predictions (`caid3/`), the metapredict 3.1.0 predictions (`metapredict_caid/`) and the PDB-disorder-monomer labels (`pdb_cache/monomer_labels.json`) from the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Each panel notebook writes `figures/figure_S09<panel>.pdf`; `figure_S09.ipynb` draws the same four panels at the same size into `figures/figure_S09.pdf`.

**The benchmark:** ordered = resolved with electron density in a PDB structure of the isolated, monomeric protein; disordered = inside such a construct but with no electron density; everything else excluded. 32,904 ordered and 3,807 disordered residues over 134 targets (8.6 : 1). Metrics are pooled over every labelled residue, as for Figure S3. The labels are built by `../../data/pdb_cache_pipeline/run_pdb_coverage.py` from PDBe (see `../../data/README.md`) and use no DisProt disorder annotations.

**Result with metapredict 3.1.0** (rank of 66): V3 AUC 0.791 (11), APS 0.299 (12), Fmax 0.412 (9), MCC 0.340 (10); V2 0.759 (19), 0.255 (21), 0.359 (18), 0.275 (20); V1 0.693 (43), 0.206 (42), 0.284 (46), 0.185 (46).

**Differences from the original figure:** metapredict from 3.1.0; one panel per metric; no suptitle.
