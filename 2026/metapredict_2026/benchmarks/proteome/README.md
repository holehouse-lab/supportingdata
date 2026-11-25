# Proteome benchmark (Figures 1E, S5)

`benchmark_metapredict.py` times one metapredict network (`--version 1|2|3`) on one device (`--device cpu|mps|gpu|auto`) predicting the human proteome in a single call, and writes a one-row TSV, `proteome_<V>_<device>_scores_<host>_<UTC>.tsv`, to `../../data/benchmarks/proteome/`.

- **Timing:** wall-clock time (`time.perf_counter`) of the prediction call alone. Reading the FASTA file and an untimed warm-up on 4 sequences are excluded, and the GPU is synchronised before the clock stops. Throughput is residues divided by this time.
- **Batch size:** not set, so each run uses metapredict 3.1's default for the device (256 on CPU and CUDA, 512 on MPS). The value used is recorded.
- **Memory:** a background thread samples the process's resident set size every 0.1 s. Peak GPU memory is `torch.cuda.max_memory_allocated` on CUDA and the largest sampled `torch.mps.driver_allocated_memory` on MPS. Memory columns are in GiB.
- **`--device auto`:** times the device metapredict picks itself (CUDA first; otherwise CPU for V1 and V2, MPS for V3).
- **`--limit N`:** predicts only the first N sequences, for quick checks. The figures drop limited runs when full runs exist.
