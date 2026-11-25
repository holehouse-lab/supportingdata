# Figure 1 — metapredict V3 is fast and accurate

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_1A.ipynb` | CAID1 Disprot-PDB: predictors ranked by MCC. |
| B | `figure_1B.ipynb` | CAID2 Disorder-PDB: predictors ranked by MCC (39 predictors). |
| C | `figure_1C.ipynb` | CAID3 Disorder-PDB: top 40 of 66 predictors ranked by MCC. |
| D | `figure_1D.ipynb` | Msn2 (P33748) disorder profiles from metapredict V1, V2 and V3, with the DNA-binding domain shaded. |
| E | `figure_1E.ipynb` | Time to predict the human proteome (20,652 sequences) with V1/V2/V3 on each device. |

`caid_mcc.py` reads the CAID files, computes MCC and draws A–C; `perf_plots.py` loads and plots the benchmark TSVs for E; `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. Each writes `figures/figure_1<panel>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Inputs, all in `../../data/` (see its README):

- A–C: `caid_references/` (reference sets; C also takes its list of disorder predictors from `caid3_disorder_pdb_results.csv`), `caid2/` and `caid3/` (CAID predictions) and `metapredict_caid/` (metapredict 3.1.0 predictions).
- E: `benchmarks/proteome/`, written by `../../benchmarks/proteome/benchmark_metapredict.py`.

**MCC scoring (A–C, as in the Methods)**

- metapredict: a residue is disordered if it lies inside a DisorderedDomain (metapredict 3.1.0 defaults). MCC for V1/V2/V3 is 0.623/0.715/0.714 on CAID1, 0.692/0.768/0.769 on CAID2 and 0.688/0.741/0.756 on CAID3 (V3 ranks 7th of 66).
- A: the other predictors' values are the published CAID1 MCCs (OpenEBench, https://openebench.bsc.es/benchmarking/OEBC011/OEBX0110000002/participants, column `dataset_default_thr_mcc`), written in the notebook.
- B: predictors are scored on their submitted binary calls. Score-only predictors are included only where a published cutoff exists (AlphaFold-rsa 0.581, AlphaFold-disorder 0.312, RONN 0.5).
- C: predictors are scored on their submitted binary calls where available. The 19 score-only predictors are cut at 0.5, except the AlphaFold2/AlphaFold3 RSA and pLDDT tracks (0.581 and 0.312).

**Timing (E):** M3 Max CPU 4.6/10/36 s for V1/V2/V3, Xeon Silver 4210R CPU 15/38/122 s, M3 Max MPS 10/20/21 s, RTX A4000 3.0/4.4/4.3 s.

**Notebook output vs the assembled figure.** A–D are used as drawn. For E, the notebook draws vertical bars labelled by machine; in the assembled figure the panel was re-laid in Illustrator with horizontal bars and retyped hardware labels. The values are the same.
