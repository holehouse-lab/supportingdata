# Batch vs sequential (Figure S7)

`run_batch_vs_sequential.py` draws random subsets of the human proteome (10 to 15,000 sequences, seeded) and times each subset two ways for each network:

- **sequential:** one `predict_disorder(sequence)` call per sequence. metapredict always predicts a single string on the CPU.
- **batch:** one `predict_disorder(list_of_sequences)` call on `--device` (default `auto`, i.e. metapredict's own choice for that network) at the default batch size.

It writes `../../data/benchmarks/batch_vs_sequential/batch_vs_sequential_<host>_<UTC>.tsv`, rewritten after each network so an interrupted run keeps its data. `--replicates` draws more independent subsets per size.

Converted from the 2025 notebook `figures/supplemental_figure_NEW_batch_vs_individual/`. The subset sizes and the two timing functions are unchanged; it is now a script that loops over V1–V3, uses a seed, synchronises the GPU, and writes TSV.
