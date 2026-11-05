# Figure 1 — metapredict V3 accuracy and speed

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_1A.ipynb` | CAID1 Disprot-PDB: predictors ranked by MCC. |
| B | `figure_1B.ipynb` | CAID2 Disorder-PDB: predictors ranked by MCC. |
| C | `figure_1C.ipynb` | CAID3 Disorder-PDB: top 40 predictors ranked by MCC. |
| D | `figure_1D.ipynb` | Msn2 (P33748) disorder profiles from metapredict V1, V2 and V3. |
| E | `figure_1E.ipynb` | Time to predict the human proteome for V1/V2/V3 across hardware. |

`caid_mcc.py` holds the CAID reading, MCC and ranking-panel functions for A–C, `perf_plots.py` loads and plots the benchmark TSVs for E, and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. A–C read the CAID reference sets (`caid_references/`), the CAID prediction files and the metapredict 3.1.0 CAID predictions from the shared data directory (`../../data`, see its README); 1C also takes its list of disorder predictors (names only) from `caid_references/caid3_disorder_pdb_results.csv`. E reads `../../data/benchmarks/proteome/`, written by `../../benchmarks/proteome/benchmark_metapredict.py` (see `../../benchmarks/README.md`). Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Each notebook writes `figures/figure_1<panel>.pdf`.

**Differences from the published panel**
- metapredict is scored from metapredict 3.1.0 predictions with 3.1.0's default binary calls (residues inside predicted IDRs are disordered), not with a 0.42/0.5 per-residue cutoff. MCC: CAID1 V1/V2/V3 = 0.623/0.715/0.714; CAID2 0.692/0.768/0.769; CAID3 0.688/0.741/0.756.
- C (CAID3) is new, D now includes V2, and E is new data from metapredict 3.1.0.
- Bars are solid black and the dataset name is written in the axes rather than as a title; A–C share one MCC axis. DFLpred (CAID1, MCC −0.029) sits below zero and has no visible bar, as before.
- `n` counts every bar shown, metapredict included.

**Values hard-coded in the notebooks.** Figure 1A uses the published CAID1 MCC values for the competing predictors (OpenEBench, https://openebench.bsc.es/benchmarking/OEBC011/OEBX0110000002/participants, column `dataset_default_thr_mcc`), because CAID1 per-residue predictions are not distributed in a usable form. The metapredict values in 1A are computed from 3.1.0 predictions.
