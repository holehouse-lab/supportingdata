# Figure S6 — prediction time vs batch size

Time to predict the human proteome (20,652 sequences) at batch sizes of 32 to 1,024, one line per machine and device (CPU black, MPS red, CUDA blue). Open squares mark the metapredict 3.1.0 default batch size for that device.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S06A.ipynb` | V1 (carries the legend). |
| B | `figure_S06B.ipynb` | V2. |
| C | `figure_S06C.ipynb` | V3. |
| A–C | `figure_S06.ipynb` | The three panels side by side with panel letters (`figures/figure_S06.pdf`). |

`perf_plots.py` loads and plots the benchmark TSVs and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. The notebooks read every TSV in `../../data/benchmarks/batch_size/`, written by `../../benchmarks/batch_size/run_batch_size_sweep.py` (see `../../benchmarks/README.md`).

**Notes.** CPU points are the median of three sweeps; GPU points are single runs. V3 ran out of the RTX A4000's 16 GB at batch sizes of 512 and 1,024, so its CUDA line stops at 256. Scores at each batch size differ from those at batch size 32 by at most 6.1 × 10⁻⁶ (`max_abs_diff_vs_min_batch` in the TSVs).
