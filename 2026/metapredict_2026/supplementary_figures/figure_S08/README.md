# Figure S8 — streaming prediction keeps memory bounded

`predict_disorder` holds the whole dataset and every prediction in memory, so its peak RAM grows with the dataset. `predict_disorder_stream` reads, predicts and returns 20,000 sequences at a time, so its peak RAM stays roughly flat. Both are run with metapredict V3 at batch size 32 on synthetic datasets of 20,000 to 2,880,000 random sequences (7.0 million to 1.0 billion residues), on an M3 Max CPU, a Xeon Silver 4210R CPU and an RTX A4000 GPU.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S08A.ipynb` | Peak RAM vs dataset size, with the `predict_disorder` trend on the M3 Max extrapolated to 16 and 64 GB. |
| B | `figure_S08B.ipynb` | Throughput vs dataset size. |
| A–B | `figure_S08.ipynb` | Both panels side by side with panel letters (`figures/figure_S08.pdf`). |

`perf_plots.py` loads and plots the benchmark TSVs and `figtools.py` is the shared style.

**Running:** open a notebook from this directory and run all cells. The notebooks read one sweep per machine and device from `../../data/benchmarks/stream/`, written by `../../benchmarks/stream/benchmark_stream_memory.py` (see `../../benchmarks/README.md`). The full sweep is used where there is one, otherwise the quick sweep; the files used are printed.

**Numbers in the legend**

- Memory is in GiB (2³⁰ bytes); for the GPU it is host RAM only. `predict_disorder` reaches 7.2–7.4 GiB at 1 billion residues. The M3 Max fit is 0.41 GiB + 6.91 GiB per billion residues, which reaches 16 GiB at 2.25 billion residues and 64 GiB at 9.2 billion.
- `predict_disorder_stream` uses 0.57–0.69 GiB on the M3 Max, 0.93–1.03 GiB on the A4000 and 0.84–1.46 GiB on the Xeon CPU.
- The fraction of residues in single-length batches (16% at 20,000 sequences, 96% at 2.88 million) is printed by `../../benchmarks/stream/single_length_batch_fraction.py`.
- The test with oneDNN's cache turned off was run separately on the Xeon and is not part of this repository.
