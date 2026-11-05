# Streaming memory benchmark (Figure S8)

`benchmark_stream_memory.py` compares the peak RAM and throughput of `metapredict.predict_disorder` (holds the whole dataset and every prediction at once) with `metapredict.predict_disorder_stream` (reads lazily and holds one chunk) on synthetic datasets of increasing size.

- **Isolation:** each measurement runs in a fresh subprocess (`stream_memory_worker.py`, copied unchanged), so peak RSS cannot leak between points.
- **Fixed settings:** both functions run V3 at batch size 32 on the chosen device (`--device`, CPU by default; Figure S8 also shows a CUDA sweep), so the activation-memory floor is small and identical, and only the dataset size varies. The chunk size is 20,000, the metapredict 3.1 default (the original sweep used 5,000).
- **Inputs:** synthetic sequences have log-normal lengths (median ~300 aa, clamped to 40–2,500) from a fixed seed, so the FASTA files in `../../data/benchmarks/stream_synthetic_data/` (not tracked) are byte-reproducible and rebuilt when missing.

Rows are appended to `../../data/benchmarks/stream/stream_<host>_<version>_<device>_bs<batch>_chunk<chunk>[_quick].tsv`. Re-running the same configuration resumes from that file.

- `python benchmark_stream_memory.py --quick` runs 3 datasets (20k–200k sequences; ~15 min on an M3 Max, ~35 min on a Xeon Silver 4210R CPU).
- `python benchmark_stream_memory.py` runs the full sweep of 7 datasets up to ~1 billion residues (~4.4 h on an M3 Max, ~10.7 h on a Xeon Silver 4210R CPU, ~1.7 h on an RTX A4000 with `--device cuda`; ~2.2 GB of synthetic FASTA).

Adapted from `figures/stream_performance/`. The method is unchanged; the output is now TSV with metapredict provenance, and memory is stored in GiB.
