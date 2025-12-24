# About

This directory has three subdirectories. Each of which has STARLING ensembles for the same set of sequences defined in the `../all_comparison_seqs.fasta` file.

* `cuda_rep1/` - GPU-generated ensembles (one replicate)
* `cuda_rep2/` - GPU-generated ensembles (another replicate)
* `mps/` - MPS-generated ensembles

These are basically identical, but we check this by comparing ensembles generated with GPUs or MPS (or independent replicas).