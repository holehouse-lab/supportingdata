# metapredict performance benchmarks

Scripts that measure how fast metapredict predicts disorder and how much memory it uses. Each script writes TSV files to its own subdirectory of `../data/benchmarks/` (the repository's data directory), and the figure notebooks read every TSV they find there. You can therefore run the scripts on several machines, copy the TSVs back into those directories, and re-run the notebooks; new machines appear in the figures without any code changes.

All figures use **metapredict 3.1.x**. Every row records `metapredict_version`, `metapredict_path` (the install that was imported) and `metapredict_git_commit`, and the notebooks warn about and drop rows from any other version.

| Directory | Script | Figure | What it measures |
| --- | --- | --- | --- |
| `proteome/` | `benchmark_metapredict.py` | 1E, S5 | Time, throughput, peak RAM and GPU memory to predict the human proteome, one version × device per run |
| `batch_size/` | `run_batch_size_sweep.py` | S6 | Prediction time vs batch size (32–1024), per version and device |
| `batch_vs_sequential/` | `run_batch_vs_sequential.py` | S7 | Time vs number of sequences, batched vs one-at-a-time |
| `stream/` | `benchmark_stream_memory.py` | S8 | Peak RAM and throughput of `predict_disorder` vs `predict_disorder_stream` as datasets grow |

`../data/benchmarks/UP000005640_9606.fasta` is the human reference proteome (UniProt UP000005640, one sequence per gene; 20,652 sequences, 11,453,425 residues), used by the first three benchmarks (see `../data/benchmarks/README.md`). `bench_env.py` holds the shared provenance and TSV helpers. The stream benchmark uses synthetic sequences generated from a fixed seed.

## Running

Install metapredict 3.1 plus `psutil` and `protfasta` (`py-cpuinfo` is optional and gives cleaner CPU names), then run from this directory. The scripts read and write `../data/benchmarks/`, so put the whole repository (not just `benchmarks/`) on each benchmark machine. Close other heavy programs first, because the timings and peak-RAM numbers are sensitive to background load. `--machine-label` sets the name a machine gets in the figures; without it, the CPU name is used (CUDA GPUs are always named by the GPU). The figures group rows by this label, so give every machine its own label (e.g. the CPU model, `--machine-label "Xeon Silver 4210R"`) rather than the generic `"x86-64"` used in the examples below.

**Apple-silicon Mac (CPU + MPS):** `bash run_local_mac.sh` runs everything below in sequence (about 1–2 h). `bash run_local_mac.sh --with-full-stream` also runs the full streaming sweep, which takes about 4 more hours and is resumable. The full sweep writes ~2 GB of synthetic FASTA to `../data/benchmarks/stream_synthetic_data/`; set `METAPREDICT_STREAM_DATA_DIR=/some/local/scratch` to put it elsewhere (e.g. outside a Dropbox-synced folder).

**Linux x86-64, CPU only:**

```bash
for v in 1 2 3; do python proteome/benchmark_metapredict.py --version $v --device cpu --machine-label "x86-64"; done
python batch_size/run_batch_size_sweep.py --devices cpu --machine-label "x86-64"
python batch_vs_sequential/run_batch_vs_sequential.py --device cpu --machine-label "x86-64"
python stream/benchmark_stream_memory.py --quick --machine-label "x86-64"   # drop --quick for the full sweep
```

**NVIDIA CUDA machine:** run the CPU commands above (so the GPU speed-up in S5D can be computed), then:

```bash
for v in 1 2 3; do python proteome/benchmark_metapredict.py --version $v --device gpu --machine-label "x86-64"; done
python batch_size/run_batch_size_sweep.py --devices cuda --machine-label "x86-64"
python batch_vs_sequential/run_batch_vs_sequential.py --device cuda --machine-label "x86-64"
python stream/benchmark_stream_memory.py --device cuda --machine-label "x86-64"   # full sweep, ~1.7 h on an RTX A4000
```

**Copying results back:** copy the TSVs in each `data/benchmarks/<benchmark>/` on the remote machine into the matching directory here (e.g. `rsync -av remote:metapredict_2026/data/benchmarks/proteome/ ../data/benchmarks/proteome/`), then re-run the notebooks. File names include the hostname and a UTC timestamp, so files from different machines never collide. If a configuration has repeat runs, the notebooks plot the median.

Files whose names start with `SMOKETEST_` came from tiny smoke tests. The notebooks only plot them when no real results exist.
