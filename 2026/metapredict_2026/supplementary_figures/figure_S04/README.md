# Figure S4 — metapredict V1, V2 and V3 against the top CAID3 predictors on five proteins

| Panel | Notebook | Contents |
| --- | --- | --- |
| — | `figure_S04.ipynb` | Per-residue disorder profiles for DP04023, DP04296, DP04242, DP04253 and DP04340, each above a bar of its Disorder-PDB labels. |

`caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data and draws the panels, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells (~15 s each). Reads the CAID3 reference set (`caid_references/`), the CAID3 predictions (`caid3/`) and the metapredict 3.1.0 predictions (`metapredict_caid/`) from the shared data directory (`../../data`, see its README). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. The notebook writes `figures/figure_S04.pdf`.

**The five proteins** (150–300 residues, metapredict V3 per-protein AUC ≥ 0.98 in the original analysis, ≥ 0.97 with metapredict 3.1.0, 21–67% disordered) were chosen in the original CAID3 re-analysis to show a range of disorder content and IDR architectures. The top five CAID3 predictors are ranked by pooled AUC over the whole reference set.

**Differences from the original figure** (`fig_si_version_profiles` in the CAID3 re-analysis)
- metapredict V1–V3 are predicted with metapredict 3.1.0.
- No suptitle; the DisProt ID heads each protein and one "Disorder score" label serves all five; the 0.5 threshold is a thin black dashed line; top CAID3 methods use the fixed manuscript colours.
- The Disorder-PDB labels are a thin bar under each trace (orange disordered, green ordered, hatched grey unscored), drawn by the same `caid3_panels.trace_panel` as Figure 2C, instead of shading. DP04023's 9 unscored residues are now visible (they used to look ordered).
