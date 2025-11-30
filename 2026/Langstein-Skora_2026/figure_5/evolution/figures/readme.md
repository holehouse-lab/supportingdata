# About

This directory contains a set of figures associated with the evolution under chemical specificity.


#### `signal_noise_round_samples/` 

This subdirectory contains comparison of per-residue conservation assessed with and without chemical conservation where the "sample size" parameter varies from 2 and 28 and reflects the number of "orthologous" sequences we use to calculate the per-residue conservation in both cases. i.e. if sample size is 28, we randomly select 28 sequences from the library of sequences under selection and 28 for those not under selection and calculate per-residue conservation from those two sets. 


#### `individual_intermaps/`
This subdirectory contains 350 individual intermaps for possible IDR variants that are under selection for chemical specificty. Note we distribute these zipped up as `individual_intermaps.zip` for storage/distribution purposes.


### `REAL_` figures

For some of these figures we the shared files start with the prefix `REAL_` - e.g.:

* `REAL_conservation_histo_with_distribution.pdf` 
* `REAL_idr2_vs_rad7_noise_vs_signal.pdf`
* `REAL_linear_conservation_no_distribution.pdf`
* `REAL_linear_conservation_with_distribution.pdf`


We do this because these figures are generated using a stochastic element in the code, so if we re-run the code the generated figures (lacking the `REAL_` prefix) would not be indentical to the figures used in the paper. 
