# Figure S8 — Streaming prediction bounds memory

`predict_disorder` holds the whole dataset and every prediction in memory, so its peak RAM grows with the dataset. `predict_disorder_stream` reads lazily and holds one chunk at a time, so its peak RAM stays roughly flat. Both are shown on every machine and device benchmarked (M3 Max CPU, Xeon Silver 4210R CPU, RTX A4000 GPU); for the GPU, peak RAM is host memory and GPU memory is not included.

| Notebook | Panel |
| --- | --- |
| `figure_S08A.ipynb` | Peak RAM vs dataset size, with the full-mode trend extrapolated to 16 and 64 GB |
| `figure_S08B.ipynb` | Throughput vs dataset size |
| `figure_S08.ipynb` | A and B side by side with panel letters (`figures/figure_S08.pdf`) |

**Data.** One sweep per machine and device from `../../data/benchmarks/stream/`, written by `benchmarks/stream/benchmark_stream_memory.py` (V3, batch size 32 for both functions, chunk size 20,000, the metapredict 3.1.0 default). The full sweep is used if present, otherwise the quick sweep; the chosen files are printed. The 16/64 GB projection is fitted to the M3 Max CPU full-mode points. Earlier sweeps at chunk size 5,000 are kept in `../../data/benchmarks/stream/chunk5000/` and not plotted. `perf_plots.py` holds the helpers and `figtools.py` the shared style. PDFs are written to `figures/`.

**Run.** Open each notebook and run all cells from this directory.

**Differences from the July 2026 draft.** The figure is regenerated with metapredict 3.1.0. The draft's 3.0.2 measurements are not reused. Colours follow the manuscript palette (black/red for the two functions) rather than the V2/V3 version colours.
