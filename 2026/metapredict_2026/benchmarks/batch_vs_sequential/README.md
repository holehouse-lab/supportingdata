# Batch vs sequential (Figure S7)

`run_batch_vs_sequential.py` draws random subsets of the human proteome (10 to 15,000 sequences; seeded) and times each subset two ways for each network:

- **sequential:** one `predict_disorder(sequence)` call per sequence. metapredict always predicts a single sequence on the CPU.
- **batch:** one `predict_disorder(list_of_sequences)` call at the default batch size, on `--device` (default `auto`: metapredict's own choice for that network).

It writes `../../data/benchmarks/batch_vs_sequential/batch_vs_sequential_<host>_<UTC>.tsv`, rewritten after each network so an interrupted run keeps its data. `--replicates` draws more subsets per size.

One random generator is shared across the three networks, so V1, V2 and V3 are timed on different subsets, but with the same seed each network gets the same subsets on every machine.
