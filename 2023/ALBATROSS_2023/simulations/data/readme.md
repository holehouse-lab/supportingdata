# Sequence data
##### Last updated 2023-11-22

## About
The files in the compressed archive `all_sequences.tgz` come from the paper:

Lotthammer, J. M.<sup>\*</sup>, Ginell, G. M.<sup>\*</sup>, Griffith, D.<sup>\*</sup>, Emenecker, R. J. & Holehouse, A. S. Direct Prediction of Intrinsically Disordered Protein Conformational Properties From Sequence. 
**Nature Methods** (in press), (2023).
*<sup>\*</sup> - these authors contributed equally.*

[Preprint available here](http://dx.doi.org/10.1101/2023.05.08.539824)

If you use these sequences and/or the sequence-prediction pairings, please cite our paper and ensure this data (via Zenodo) is directly referenced ([doi:10.5281/zenodo.10198621](https://zenodo.org/records/10198621)).

## Structure
The `all_sequences.tgz` file is a compressed directory with the following structure:

	.
	├── test
	│   ├── asph_nat_meth_test.tsv
	│   ├── prefactor_nat_meth_test.tsv
	│   ├── re_nat_meth_test.tsv
	│   ├── rg_nat_meth_test.tsv
	│   ├── scaled_re_nat_meth_test.tsv
	│   ├── scaled_rg_nat_meth_test.tsv
	│   └── scaling_exp_nat_meth_test.tsv
	└── training
	    ├── asph_bio_synth_training_data_cleaned_05_09_2023.tsv
	    ├── prefactor_bio_synth_training_data_cleaned_05_09_2023.tsv
	    ├── re_bio_synth_training_data_cleaned_05_09_2023.tsv
	    ├── rg_bio_synth_training_data_cleaned_05_09_2023.tsv
	    ├── scaled_re_bio_synth_training_data_cleaned_05_09_2023.tsv
	    ├── scaled_rg_bio_synth_training_data_cleaned_05_09_2023.tsv
	    └── scaling_exp_bio_synth_training_data_cleaned_05_09_2023.tsv
	    
## Data sets
There are two subdirectories; **test** and **training**. The sequences in each of these were generated as described in the methods. To re-iterate for convenience:

#### Training data
Using the IDR design package [GOOSE](https://github.com/idptools/goose), we assembled a library of chemically diverse synthetic disordered proteins [1]. Sequences varied charge, hydropathy, charge patterning, and amino acid composition and were between 10 and 750 residues long. We generated 23,127 disordered protein sequences across a diverse sequence space (see Supplementary Information). In addition to the synthetic sequence library, we curated a set of 19,075 naturally occurring IDRs by randomly sampling disordered proteins ranging in length from 10-750 residues from one of each of the following proteomes: *Homo sapien*, *Mus musculus*, *﻿Dictyostelium discoideum*, *Escherichia coli*, *Drosophilia melanogaster*, *Saccharomyces cerevisiae*, *Neurospora crassa*, *Schizosaccharomyces pombe*, *Xenopus laevis*, *Caenorhabditis elegans*, *Arabidopsis thaliana*, and *Danio rerio*. All annotated IDRs from the proteomes mentioned above are available at [https://github.com/holehouse-lab/shephard-data/tree/main/data/proteomes](https://github.com/holehouse-lab/shephard-data/tree/main/data/proteomes).

#### Test data
To prepare a test set to accurately assess the true generalization error for each of our ALBATROSS predictors, we randomly selected an additional set of 2,501 biological IDRs from one of the aforementioned proteomes. To ensure that any newly selected biological sequences were distinct from those seen during training, we applied [CD-HIT](https://sites.google.com/view/cd-hit) with default parameters to remove sequences with >20% similarity, leaving 2,306 biological IDRs in our test set. We also designed and simulated 3,731 synthetic disordered protein sequences using the same design parameters as in our training set. All sequences generated were between 10 and 750 residues in length. The ALBATROSS test set we used to assess model accuracy consisted of 6,037 disordered protein sequences.
	    
## Data types	    
There are seven types of analysis for each sequence in the training and test datasets.

These are:

1. `rg_` - **Radius of gyration**: The overall ensemble-average radius of gyration of the ensemble.
2. `re_` - **End-to-end distance**: The overall ensemble-average end-to-end distance associated with the ensemble.
3. `scaled_rg` - **Radius of gyration** divided by the square root of the length: Same as the radius of gyration except divided by sqrt(N). This division is a normalization factor in polymer space (because polymer dimensions scale as $N^{\nu}$ and $N^{0.5} = \sqrt N$. This basically removes length as a determinant of chain dimensions meaning we effectively regularize the radius of gyration according to polymer theory. For more information, recommended *Polymer Physics* by Rubinstein and Colby (both here, but in general, for all life events).
4. `scaled_re` - **End-to-end distance** divided by the square root of the length. Same logic as the `scaled_Rg` data.
5.  `asph_` - **Asphericity**: The overall ensemble asphericity, as defined and calculated by [SOURSOP](https://soursop.readthedocs.io/en/latest/). For mathematical definition see equation 1 in [Lalmansingh et al. 2023](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00190).
6. `prefactor_` - Prefactor ($A$) in fitting the internal scaling profiles to an exponential function of the form $R_{ij} = A|i-j|^{\nu}$, where $R_{ij}$ reflects all possible intramolecular distances, A is the prefactor in question, |i-j| reflects all possible sequence separations where i and j here match those in $R_{ij}$ and $\nu$ is the apparent scaling exponent obtained for this fitting. For a more detailed discussion on the process of extracting an apparent scaling exponent from ensembles of finite-sized chains see [Peran & Holehouse et al. 2019](https://www.pnas.org/doi/suppl/10.1073/pnas.1818206116). Basically, if the chain is well-described by a homopolymer then $R_e = AN^{\nu}$.
7. `scaling_exp_` - Scaling exponent ($\nu$) in fitting the internal scaling profiles to an exponential function of the form $R_{ij} = A|i-j|^{\nu}$, where $R_{ij}$ reflects all possible intramolecular distances, A is the prefactor in question, |i-j| reflects all possible sequence separations where i and j here match those in $R_{ij}$ and $\nu$ is the apparent scaling exponent obtained for this fitting. For a more detailed discussion on the process of extracting an apparent scaling exponent from ensembles of finite-sized chains see [Peran & Holehouse et al. 2019](https://www.pnas.org/doi/suppl/10.1073/pnas.1818206116). Basically, if the chain is well-described by a homopolymer then $R_e = AN^{\nu}$. This scaling exponent can be used to estimate the solvent quality
	    
	    
## Data file structure	    
	    
Each file is a tab-separated file with three columns:

1. **Sequence ID**; a unique identifier for each sequence
2. **Sequence**; amino acid sequence
3. **Value**; the specific numerical value associated with the analysis relevant to that file as a result of Mpipi-GG simulations, described above.

To reiterate, the "value" column here is the actual simulated value. 

