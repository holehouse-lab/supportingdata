# metapredict performance benchmarks

Scripts that measure how fast metapredict predicts disorder and how much memory it uses (Figures 1E and S5–S8). Each script writes TSVs to its own subdirectory of `../data/benchmarks/`, and the figure notebooks plot every TSV they find there, so results from a new machine appear in the figures without code changes.

| Directory | Script | Figure | What it measures |
| --- | --- | --- | --- |
| `proteome/` | `benchmark_metapredict.py` | 1E, S5 | Time, throughput, peak RAM and GPU memory to predict the human proteome, one network and device per run. |
| `batch_size/` | `run_batch_size_sweep.py` | S6 | Prediction time at batch sizes of 32–1,024, per network and device. |
| `batch_vs_sequential/` | `run_batch_vs_sequential.py` | S7 | Time vs. number of sequences, batched vs. one per call. |
| `stream/` | `benchmark_stream_memory.py` | S8 | Peak RAM and throughput of `predict_disorder` vs `predict_disorder_stream` as datasets grow. |

`bench_env.py` holds the shared paths, provenance, and TSV helpers. The first three benchmarks predict the human reference proteome (`../data/benchmarks/UP000005640_9606.fasta`; 20,652 sequences, 11,453,425 residues); the stream benchmark generates synthetic sequences from a fixed seed. Every row records the hardware, operating system, Python, PyTorch, and metapredict version (with path and git commit); the notebooks drop rows from any metapredict other than 3.1.x.

## Machines in the paper

| Label | Hardware | Software |
| --- | --- | --- |
| M3 Max | 14-inch MacBook Pro (2023), Apple M3 Max, 128 GB; GPU runs use MPS | macOS 26.5.1, Python 3.12.13, PyTorch 2.10.0 |
| Xeon Silver 4210R / RTX A4000 | Linux workstation, Intel Xeon Silver 4210R (10 cores), 64 GB, NVIDIA RTX A4000 (16 GB); GPU runs use CUDA | Ubuntu 20.04, Python 3.14.8, PyTorch 2.14.1 (CUDA 13.0) |

## Running

Install metapredict 3.1 plus `psutil` and `protfasta` (`py-cpuinfo` is optional and gives cleaner CPU names), and run from this directory. The scripts read and write `../data/benchmarks/`, so each benchmark machine needs this directory and `data/benchmarks/`. Close other heavy programs first: the timings and peak-RAM numbers are sensitive to background load. `--machine-label` sets the name a machine gets in the figures (default: the CPU name; CUDA GPUs are always named by the GPU). The figures group rows by this label, so give every machine its own.

**Apple-silicon Mac (CPU + MPS):** `bash run_local_mac.sh` runs everything below in sequence (1–2 h). `bash run_local_mac.sh --with-full-stream` adds the full streaming sweep (about 4 more hours, resumable).

**Linux x86-64, CPU:**

```bash
for v in 1 2 3; do python proteome/benchmark_metapredict.py --version $v --device cpu --machine-label "x86-64"; done
python batch_size/run_batch_size_sweep.py --devices cpu --machine-label "x86-64"
python batch_vs_sequential/run_batch_vs_sequential.py --device cpu --machine-label "x86-64"
python stream/benchmark_stream_memory.py --quick --machine-label "x86-64"   # drop --quick for the full sweep
```

**NVIDIA CUDA machine:** run the CPU commands above (so the GPU speed-up in Figure S5D can be computed), then:

```bash
for v in 1 2 3; do python proteome/benchmark_metapredict.py --version $v --device gpu --machine-label "x86-64"; done
python batch_size/run_batch_size_sweep.py --devices cuda --machine-label "x86-64"
python batch_vs_sequential/run_batch_vs_sequential.py --device cuda --machine-label "x86-64"
python stream/benchmark_stream_memory.py --device cuda --machine-label "x86-64"   # full sweep, ~1.7 h on an RTX A4000
```

**Copying results back:** copy the TSVs from each `data/benchmarks/<benchmark>/` on the remote machine into the same directory here, then re-run the notebooks. File names include the hostname and a UTC timestamp, so files from different machines never collide. Repeat runs of one configuration are plotted as their median. Files starting with `SMOKETEST_` come from tiny test runs and are plotted only when no real results exist.

The full stream sweep writes ~2 GB of synthetic FASTA to `../data/benchmarks/stream_synthetic_data/`; set `METAPREDICT_STREAM_DATA_DIR` to put it elsewhere.
