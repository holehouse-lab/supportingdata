# Figure 2 — metapredict V3 against the CAID3 field

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_2A.ipynb` | ROC curves: metapredict V1/V2/V3 and the five CAID3 predictors with the highest AUC; the other 58 CAID3 disorder predictors in grey. |
| B | `figure_2B.ipynb` | Precision–recall curves, as in A. |
| C | `figure_2C.ipynb` | Per-residue traces of metapredict V3 and the top five CAID3 predictors, each above its Disorder-PDB bar. Draws four proteins; the published panel shows three (DP04242, DP04068, DP03933). |
| D (left) | `figure_2D_left.ipynb` | The five proteins with the highest metapredict V3 per-protein AUC, with PDB structural evidence. The published panel shows three (DP04263, DP04172, DP04297). |
| D (right) | `figure_2D_right.ipynb` | The five proteins with the lowest per-protein AUC, drawn the same way. The published panel shows three (DP04228, DP04249, DP04214). |

**Notebook output vs the assembled figure.** A and B are used as drawn. The notebooks draw more proteins than the figure shows: four in C (three shown) and five in each half of D (three shown). The others were removed in Illustrator, so `figures/figure_2C.pdf`, `figure_2D_left.pdf` and `figure_2D_right.pdf` are longer than the published panels. Set `TRACE_TARGETS` (C) or `N_EXAMPLES` (D) in the notebooks to change what is drawn.

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code, `caid3_panels.py` loads the data and draws A–C, `ground_truth_panels.py` draws the four-track blocks in D, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (15–30 s each). Each writes `figures/<notebook name>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid_references/caid3_disorder_pdb.fasta`, `caid3/`, `metapredict_caid/` and, for D, `pdb_cache/` from `../../data/` (see its README).

**Scoring (as in the Methods)**

- The 63 CAID3 disorder predictors (binding and linker predictors, listed in `predictor_classes.py`, are excluded) plus metapredict V1–V3 are scored on CAID3 Disorder-PDB: 319 proteins, 31,401 disordered and 67,838 ordered residues, pooled. Each predictor is scored on the proteins it predicted.
- AUC and APS come from scikit-learn, Fmax is the highest F1 over all cutoffs, and MCC is taken at the Fmax cutoff. The notebooks print each highlighted predictor's value and rank.
- Per-protein AUC (D) is the ROC AUC of the V3 score over one protein's scored residues. Only proteins of 40–500 residues with at least 15 disordered and 15 ordered scored residues are eligible (116 proteins; 27 have an AUC of 1.00 and are ranked by their number of disordered residues).
- PDB track (D): each residue is labelled as resolved in a monomer, a homo-oligomer or only a hetero-complex, in a construct but with missing density, or in no PDB construct. The PDB codes behind each block are printed by the notebook.
- IDR shading (D) uses metapredict's default domain calling (0.5 threshold, minimum IDR 12 residues, minimum folded domain 50, gap closure 10).

**Result:** metapredict V3 has AUC 0.929 (14th of 66), APS 0.902 (7th), Fmax 0.828 (8th) and MCC 0.761 (7th); see Figure S3 for V1 and V2.
