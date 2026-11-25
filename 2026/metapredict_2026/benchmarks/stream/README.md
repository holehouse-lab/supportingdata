# Streaming memory benchmark (Figure S8)

`benchmark_stream_memory.py` compares the peak RAM and throughput of `metapredict.predict_disorder` (holds the whole dataset and every prediction) with `metapredict.predict_disorder_stream` (reads, predicts and returns one chunk at a time) on synthetic datasets of increasing size.

- **Datasets:** sequence lengths are log-normal (μ = 5.7, σ = 0.55 on the natural-log scale; median ~300 residues), rounded down and limited to 40–2,500 residues; residues are drawn uniformly from the 20 amino acids. The full sweep has 20,000, 80,000, 250,000, 500,000, 720,000, 1,440,000 and 2,880,000 sequences (7.0 million to 999 million residues). The generator is re-seeded for each dataset, so the files are identical on every machine.
- **Settings:** V3 at batch size 32 for both functions, on `--device` (CPU by default), with a chunk size of 20,000 (the metapredict 3.1 default) for streaming.
- **Measurement:** each point runs in a fresh subprocess (`stream_memory_worker.py`), which samples its resident set size every 4 ms. Memory is stored in GiB.

Rows are appended to `../../data/benchmarks/stream/stream_<host>_<version>_<device>_bs<batch>_chunk<chunk>[_quick].tsv`; re-running a configuration resumes from that file. The FASTA files go in `../../data/benchmarks/stream_synthetic_data/` (~2 GB for the full sweep) and are rebuilt when missing.

- `python benchmark_stream_memory.py --quick`: 3 datasets of 20,000 to 200,000 sequences (~15 min on an M3 Max, ~35 min on a Xeon Silver 4210R CPU).
- `python benchmark_stream_memory.py`: the full sweep (~4.4 h on an M3 Max, ~10.7 h on a Xeon Silver 4210R CPU, ~1.7 h on an RTX A4000 with `--device cuda`).
- `python single_length_batch_fraction.py`: prints, for each dataset, the fraction of residues that fall in batches where every sequence has the same length (16% at 20,000 sequences and 96% at 2.88 million for `predict_disorder`; ~16% at every size for `predict_disorder_stream`). It redraws the lengths from the seed and needs no prediction.

The Figure S8 legend also describes a 200,000-sequence test on the Xeon with oneDNN's primitive cache turned off; that test is not part of these scripts.
