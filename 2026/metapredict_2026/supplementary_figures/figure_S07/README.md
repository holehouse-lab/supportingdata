# Figure S7 — batch vs sequential prediction

Time to predict random subsets of 10 to 15,000 human proteins, either one sequence per call (sequential, red; always on the CPU) or in one batched call at the default batch size (black on a CPU, blue on a GPU). Dashed lines are linear fits.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S07A.ipynb` | V1. |
| B | `figure_S07B.ipynb` | V2. |
| C | `figure_S07C.ipynb` | V3. |
| A–C | `figure_S07.ipynb` | The three panels side by side with panel letters (`figures/figure_S07.pdf`). |

`perf_plots.py` loads and plots the benchmark TSVs and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. The notebooks read every TSV in `../../data/benchmarks/batch_vs_sequential/`, written by `../../benchmarks/batch_vs_sequential/run_batch_vs_sequential.py` (see `../../benchmarks/README.md`).

**Notes.** On the MacBook Pro, batched calls use each network's default device (CPU for V1 and V2, MPS for V3). The workstation was run twice, batching on the CPU and on the GPU; its sequential times are the median of the two runs.
