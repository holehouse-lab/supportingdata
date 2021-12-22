## Simulation data
###### Last updated 2021-12-22

This directory contains output from the all-atom simulations performed in Jankowski *et al.* 

Output from all-atom simulations performed using the standard ABSINTH forcefield are provided under the following directories:

* `/WT` (wildtype sequence)
* `/AA` (sequence with double arginine mutated to alanine)
* `/HH` (sequence with double arginine mutated to neutral histidine)
* `/HH_pro` (sequence with double arginine mutated to protonated histidine)

Within these directories, five files exist:

* `RG.csv` - instantaneous radius of gyration for every confomer
* `DSSP_H.csv` - per-residue mean helicity across the ensemble, as calculated by the DSSP algorithm
* `SASA_SC_mean.csv` - per-residue solvent accessible surface area (in angstrom squared) for the sidechains from each residue
* `start.pdb` - Topology file defining atomic positions 
* `traj.xtc` - Sub-sampled trajectory file 1000 confomers. If the full trajectory files are requested they will be deposited and shared on Zenodo.

For excluded volume (EV) simulations, analogous directories can be found under `/EV_simulations`.

