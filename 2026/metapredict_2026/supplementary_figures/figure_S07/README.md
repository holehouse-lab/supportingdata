# Figure S7 — Batch vs one-at-a-time prediction

Time to predict random subsets of the human proteome (10–15,000 sequences), either one sequence per call (sequential, red; always on the CPU) or in a single batched call (black on CPU, blue on GPU; default batch size, on metapredict 3.1.0's default device unless the run pinned one, e.g. the Xeon CPU series on `aviation`). Dashed lines are linear fits.

| Notebook | Panel |
| --- | --- |
| `figure_S07A.ipynb` | V1 |
| `figure_S07B.ipynb` | V2 |
| `figure_S07C.ipynb` | V3 |
| `figure_S07.ipynb` | A–C side by side with panel letters (`figures/figure_S07.pdf`) |

**Data.** Every TSV in `../../data/benchmarks/batch_vs_sequential/`, written by `benchmarks/batch_vs_sequential/run_batch_vs_sequential.py` (see `benchmarks/README.md`). `perf_plots.py` holds the helpers and `figtools.py` the shared style. PDFs are written to `figures/`.

**Run.** Open each notebook and run all cells from this directory.

**Differences from the 2025 draft.** That draft was one panel for a single version, with the timings made inside the notebook. The timing is now a seeded script covering all three versions, with one panel per version.
