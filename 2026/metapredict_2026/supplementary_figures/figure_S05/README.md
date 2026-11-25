# Figure S5 — speed and memory across versions and hardware

metapredict V1, V2 and V3 predicting the human proteome (20,652 sequences, 11,453,425 residues) at the default batch size for each device (256 on CPU and CUDA, 512 on MPS). These are the same runs as Figure 1E.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S05A.ipynb` | Throughput (residues/s, log scale). |
| B | `figure_S05B.ipynb` | Peak RAM (resident set size of the Python process). |
| C | `figure_S05C.ipynb` | Peak GPU memory. |
| D | `figure_S05D.ipynb` | GPU speed-up over the CPU of the same machine. |

`perf_plots.py` loads and plots the benchmark TSVs and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. Each writes `figures/figure_S05<panel>.pdf`. The notebooks read every TSV in `../../data/benchmarks/proteome/`, written by `../../benchmarks/proteome/benchmark_metapredict.py` (see `../../benchmarks/README.md`). Rows from metapredict versions other than 3.1.x are dropped with a warning.

**Result:** GPU speed-up is 4.9×, 8.6× and 28× for V1, V2 and V3 on the RTX A4000, and 0.4×, 0.5× and 1.7× on the M3 Max (MPS), where the CPU is faster for V1 and V2. Peak RAM for V3 on the CPU is 16.0 GiB (M3 Max) and 11.8 GiB (Xeon), set by the longest sequence (titin, 34,350 residues).
