# Batch-size sweep (Figure S6)

`run_batch_size_sweep.py` predicts the human proteome at batch sizes of 32, 64, 128, 256, 512 and 1,024 for each network (`--versions`) on each device (`--devices cpu mps cuda`), and appends one row per measurement to `../../data/benchmarks/batch_size/batch_size_<host>_<UTC>.tsv`.

- **Timing:** an untimed warm-up runs first, and the device is synchronised before each timer stops.
- **Memory:** cached GPU memory is released before each measurement. A batch size that runs out of memory is reported and gets no row. On the 16 GB RTX A4000, V3 does not fit at 512 or 1,024, because the first batch is padded to titin's length (34,350 residues).
- **Prediction check:** `max_abs_diff_vs_min_batch` is the largest difference in unrounded scores from the smallest batch size (at most 6.1 × 10⁻⁶), confirming that batch size does not change predictions.
- **Defaults:** `default_batch_size` is the batch size metapredict 3.1 would choose on that device; the figure marks it.
- **`--cooldown 600`:** rests 10 minutes before each measurement, to avoid thermal throttling.

In the paper, the CPU sweep was run three times on each machine (the figure plots the median) and each GPU sweep once.
