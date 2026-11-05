# Figure 2 — metapredict V3 against the CAID3 field, and where it is right and wrong

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_2A.ipynb` | ROC curves: metapredict V1/V2/V3 and the top 5 CAID3 predictors, with every other CAID3 predictor in grey. |
| B | `figure_2B.ipynb` | Precision–recall curves, as in A. Carries the predictor legend for the whole figure. |
| C | `figure_2C.ipynb` | Per-residue traces of the top 5 CAID3 predictors and metapredict V3 on four example proteins, each above its Disorder-PDB bar. Carries the Disorder-PDB key. |
| D | `figure_2D.ipynb` | The five proteins metapredict V3 predicts best (per-protein AUC), with PDB structural evidence. Carries the track names and the IDR-call/PDB key for D and E. |
| E | `figure_2E.ipynb` | The five proteins metapredict V3 predicts worst, as in D (no track names or legend; it sits to the right of D). |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data and draws A–C, `ground_truth_panels.py` draws the four-track blocks in D–E, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (each takes ~15–30 s). Needs the CAID3 reference set (`caid_references/caid3_disorder_pdb.fasta`), the CAID3 prediction files, the metapredict 3.1.0 CAID predictions and, for D–E, the PDB-coverage cache (`pdb_cache/`), all in the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Each notebook writes `figures/figure_2<panel>.pdf`.

**Layout.** The panels are drawn at print size: A and B 3.2 in wide each, C 6.45 in, D 3.5 in and E 2.9 in, so A|B and D|E each fill a 6.4 in row. Text is 7 pt (ticks, headers), 8 pt (axis labels) and 6.5 pt (legends), as in Figures 1, S10 and S11. Each colour is keyed once, so the legend says "colours as in B" for the other panels. The notebooks print what the panels no longer show: each highlighted predictor's AUC/APS and rank (A, B) and the PDB codes behind each resolved block (D, E).

**Result with metapredict 3.1.0** (66 CAID3 disorder predictors): V3 AUC 0.929 (14th), APS 0.902 (7th), Fmax 0.828 (8th), MCC 0.761 (7th); V2 0.915/0.875/0.804/0.721; V1 0.887/0.837/0.749/0.640.

**How the panels are scored.** All metrics are pooled over every covered, non-excluded residue of every target (as in the CAID3 assessment): AUC and APS from scikit-learn; Fmax as the maximum F1 over thresholds; MCC at the Fmax threshold. Only disorder predictors are ranked; binding and linker flavours (listed in `predictor_classes.py`) are excluded. Per-protein AUC (D–E) uses one target's scored residues; only targets of 40–500 residues with ≥15 disordered and ≥15 ordered scored residues are eligible. metapredict IDR shading in D–E applies metapredict's domain decomposition (0.5 threshold, 12-residue minimum IDR, 50-residue minimum folded domain, 10-residue gap closure) to the plotted scores.

**Differences from the original CAID3 analysis figures**
- metapredict is predicted with metapredict 3.1.0 rather than taken from the CAID3-distributed files (raw scores agree to ±0.001).
- Each top CAID3 method has one fixed colour in every panel (PUNCH2 green, PUNCH2-Light violet, AlphaFold-rsa teal, SPOT-Disorder2 black, AlphaFold3-rsa brown); previously the ground-truth figures used different colours for the same methods.
- Titles are removed and text is kept to a minimum: only the DisProt ID heads each protein in C, and the ID and per-protein AUC in D–E; the PDB codes are printed rather than drawn. Threshold lines are thin black dashes.
- A and B draw every CAID3 disorder predictor other than the highlighted eight in grey (previously only those below metapredict V3, which left out the methods ranked 6–13).
- The Disorder-PDB labels are a thin bar under each trace in C–E (orange disordered, green ordered, hatched grey unscored), as in Figures S4 and S10. C previously shaded disordered residues instead.
- The PDB track in D–E is drawn in greys (monomer, homo-oligomer, hetero-complex only, from dark to light) with missing density hatched pale orange, so it shares no colour with a predictor line; the NMR Cα RMSF overlay is off (pass `rmsf=load_rmsf()` to `make_composite` to bring it back).
- D/E show five proteins each, which makes them tall for a main-text panel; change `N_EXAMPLES` in the notebook to show fewer.
