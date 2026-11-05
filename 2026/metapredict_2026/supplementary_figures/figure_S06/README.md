# Figure S6 — Prediction time vs batch size

Time for each metapredict network to predict the human proteome as a function of batch size, one line per machine and device (CPU black, Apple MPS red, CUDA blue). Open squares mark the metapredict 3.1.0 default batch size for that device. These measurements motivated 3.1.0's per-device defaults: the small V1 and V2 networks are dominated by per-batch overhead on GPUs, so larger batches help a lot there.

| Notebook | Panel |
| --- | --- |
| `figure_S06A.ipynb` | V1 |
| `figure_S06B.ipynb` | V2 |
| `figure_S06C.ipynb` | V3 |
| `figure_S06.ipynb` | A–C side by side with panel letters (`figures/figure_S06.pdf`) |

**Data.** Every TSV in `../../data/benchmarks/batch_size/`, written by `benchmarks/batch_size/run_batch_size_sweep.py` (see `benchmarks/README.md`). Only sweeps over the largest input (the full proteome) are combined, and repeat measurements are summarised by their median. `perf_plots.py` holds the helpers and `figtools.py` the shared style. PDFs are written to `figures/`.

**Run.** Open each notebook and run all cells from this directory.
