# About
#### Last updated 2024-02-04

This directory contains information used for simulations associated with ALBATROSS training.

Firstly, the file `lammps_Mpipi_GGv23.in` contains an annotated lammps input file that defines the Mpipi-GGv23-r1 parameters. Note 23 = 2023 and r1 = revision 1; We do not expect to update these any time soon, but just explaining the nomenclature!

Secondly, the directory `data/` contains:

1. All sequences used for ALBATROSS training and validation
2. All simulation-derived parameters used for ALBATROSS training and validation

A detailed description of the file structure is provided in the associated README.

If the raw MD trajectories are of interest, we can work together to share these, but they constitute 100s of GBs, so this requires some effort on our end to coordinate. 


## Citation
Lotthammer, J. M.\*, Ginell, G. M.\*, Griffith, D.\*, Emenecker, R. J. & Holehouse, A. S. Direct prediction of intrinsically disordered protein conformational properties from sequences. Nat. Methods (2024). doi:10.1038/s41592-023-02159-5
