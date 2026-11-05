# Proteome benchmark (Figure 1E, Figure S5)

`benchmark_metapredict.py` times one metapredict network (`--version 1|2|3`) on one device (`--device cpu|mps|gpu|auto`) predicting disorder for the human proteome (`../../data/benchmarks/UP000005640_9606.fasta`), and writes a one-row TSV to `../../data/benchmarks/proteome/` named `proteome_<V>_<device>_scores_<host>_<UTC>.tsv`.

- **Timing:** only inference is timed. FASTA reading and the model/device initialisation (an untimed warmup on 4 sequences) are excluded, and the GPU is synchronised before the clock stops.
- **Batch size:** it is *not* set, so each run uses metapredict 3.1's default for that network and device (256 on CPU/CUDA, 512 on MPS). The resolved value is recorded in `batch_size`.
- **`--device auto`:** times the device metapredict picks by itself (V1/V2 prefer CPU over MPS, V3 prefers MPS; CUDA always first).
- **Memory:** a background thread samples peak RSS every 0.1 s. Peak GPU memory is `torch.cuda.max_memory_allocated` (CUDA) or the largest sampled `torch.mps.driver_allocated_memory` (MPS). Memory columns are in GiB.
- **Origin:** adapted from `figures/Figure_1E/benchmark_metapredict.py` (metapredict 3.0.2 era). The monitor and timed region are unchanged; the output is now TSV with metapredict provenance (version, path, git commit).

`--limit N` benchmarks only the first N sequences, for quick checks. Figures drop limited runs whenever full-proteome runs exist.
