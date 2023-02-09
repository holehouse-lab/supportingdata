## Additional Supporting Information
###### Last updated 2023-02-09 [thanks Mina!]

This directory contains information for the manuscript:

Holehouse, A. S., Ginell, G. M., Griffith, D., & Böke, E. (2021). Clustering of Aromatic Residues in Prion-like Domains Can Tune the Formation, State, and Organization of Biomolecular Condensates. Biochemistry, 60(47), 3566–3581. [link here](https://doi.org/10.1021/acs.biochem.1c00465). 


## Contents

### `/Bioinformatics`
Contains data used for bioinformatic analysis in the paper. At the time of writing the IWD aromatic clustering was calculated using in-development scripts, but this can now be done using [sparrow](https://github.com/idptools/sparrow).

### `/simulations`
Contains two main files for running PIMMS simulation

* `KEYFILE.kf` - The PIMMS keyfile used to run simulations. This keyfile has 300 copies of a `FSSSSFSSSSFSSSSFSSSSFSSSSFSSSSFSSSSFSSSSFSSSSFSSSSFSSSSF` chain, but could be changed to use different chain architectures as described in the paper.
* `parameter_input.prm` - the parameter file used to run the simulation. This the only values changed across the simulations are the sticker-spacer and spacer-spacer interactions, as highlighted in the parameter file.

Also contains a directory called `msd_fits/`  which contains mean-square displacement fits made for polymers inside droplets for sticker-sticker and sticker-spacer titrations of the uniformly spaced sticker sims (i.e. sims used in Figure 4). Note that in that figure, the average across three replicates is shown, but here we show the MSD fits for each individual simulation. Also, all fits are shown, but MSD fits are only plotted where alpha was ~1.0.


