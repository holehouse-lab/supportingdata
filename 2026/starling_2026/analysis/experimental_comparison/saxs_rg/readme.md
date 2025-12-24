# About
The code and data in these directories capture 


### Ensemble generation
Ensemble generation was done using:

	cd ensembles/
	starling ../all_comparison_seqs.fasta -b 300 -c 600 -r
	
>**NB:** Note you may need to adjust the batch-size `-b 300` arguement depending on available memory.	


### Directory structure

* `all_comparison_data.csv` - this CSV file contains 137 sequences, which are the exact same 137 sequences used previously by [Lotthammer, Ginell, Griffith et al. 2024](https://www.nature.com/articles/s41592-023-02159-5) to calibrate the ALBATROSS method. These 137 sequences are our starting point, and we filter these for those that are less than 384 residues in length, yielding 133 sequences. This file is read in and parsed by the `sequence_generation.ipynb` notebook.
* `sequence_generation.ipynb` - notebook that parses `all_comparison_data.csv` and generates as output three files (`starling_comparison_data.csv`, `all_comparison_seqs.fasta`, and `all_comparison_seqs_GS_versions.fasta`)
* `ensembles/` - this is where the STARLING ensembles for these comparisons are located. Note, we actually generated three sets of ensembles, one on a MacBook using MPS, and two using CUDA. We do this so we can compare reproducibility between independent runs (CUDA 1 vs. CUDA 2) and compare ensembles between MPS and CUDA. We see no difference across all these cross-comparisons, but we felt it was important to show explicitly!
* `average_rg_vals.ipynb` - main script that does all the analysis in this directory, generating the figures in `figures/` and the output data file in `data_out/`.
* `figures/` - figures comparing different models
* `data_out/` - CSV file with comparison $R_g$ values from the various models.

