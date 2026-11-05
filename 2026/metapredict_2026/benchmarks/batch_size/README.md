# Batch-size sweep (Figure S6)

`run_batch_size_sweep.py` predicts the human proteome at batch sizes 32–1024 for each network (`--versions`) on each available device (`--devices cpu mps cuda`). It appends one row per measurement to `../../data/benchmarks/batch_size/batch_size_<host>_<UTC>.tsv`.

- **Timing:** an off-clock warmup runs at the default batch size, and the device is synchronised before each timer stops.
- **Prediction check:** predictions use `round_values=False`. `max_abs_diff_vs_min_batch` is the largest score difference from the smallest batch size, which confirms that batch size changes speed and memory but not predictions (differences are at most ~6e-6, floating-point noise).
- **Memory:** cached GPU memory is released off the clock before each timed run, so every batch size starts from the same allocator state. A batch size that runs out of memory (GPU, MPS or CPU) is reported on screen and in a summary at the end and gets no row; the sweep carries on. On a 16 GB RTX A4000, V3 at batch 512 and 1024 does not fit, because every sequence in the first batch is padded to titin's length (34,350 aa).
- **Defaults:** `default_batch_size` is the batch size metapredict 3.1 would choose by itself on that device; the figure marks it.
- **Cooldown:** `--cooldown 600` rests 10 minutes before each measurement to avoid thermal throttling (the CPU numbers are otherwise somewhat noisy).
- **Origin:** adapted from `V3_manuscript/benchmark/batch_size/run_batch_size_sweep.py`. The batch size is now passed as `predict_disorder(batch_size=...)`. The original changed `parameters['batch_size']` in memory, which metapredict 3.1's per-device defaults override.
