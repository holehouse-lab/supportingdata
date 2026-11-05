# Figure S5 — metapredict performance across versions and hardware

Throughput, memory and GPU speed-up for metapredict V1, V2 and V3 predicting the human proteome, from the same runs as Figure 1E. This replaces the 2024 version of S5, which used hard-coded timings.

| Notebook | Panel |
| --- | --- |
| `figure_S05A.ipynb` | Throughput (residues/s, log scale) |
| `figure_S05B.ipynb` | Peak RAM (process RSS) |
| `figure_S05C.ipynb` | Peak GPU memory |
| `figure_S05D.ipynb` | GPU speed-up over the CPU of the same machine |

`perf_plots.py` holds the loading and plotting helpers (an identical copy is in figure_1, S6, S7 and S8), and `figtools.py` holds the shared style. PDFs are written to `figures/`.

**Data.** The notebooks read every TSV in `../../data/benchmarks/proteome/`, written by `benchmarks/proteome/benchmark_metapredict.py`; see `benchmarks/README.md` for how to run it on each machine. Rows from metapredict versions other than 3.1.x are dropped with a warning.

**Run.** Open each notebook and run all cells, using this directory as the working directory.

**Differences from the 2024 S5.** That figure showed hard-coded timings of "single" vs "batch" prediction on one Intel CPU and one NVIDIA A4500. These panels plot measured 3.1.0 runs at the default batch sizes on any machines in the results, and add memory and GPU speed-up.
