## Simulation Supporting Information
###### Last updated 2024-03-02

This directory describes content for simulation component of the paper "*The interaction between a disordered clock protein and its partner turns an hourglass into a persistent circadian oscillator*" by Jankowski et al. (in press at Nature Communications).

## Contents:

* `/data` - well structured simulation data for all simulations including example trajectories and full trajectory analysis data used in this paper.
* `/input_files` - Information on how to run CAMPARI simulations and input files to enable full reproducibility 
* `/manuscript_figures` - Jupyter notebooks and data for reproducing simulation analysis figures in the 

## Simulation analysis
Simulations were analyzed using [SOURSOP](https://soursop.readthedocs.io/) [1], all analysis are standard functions in SOURSOP so the code for analyzing the simulations is not actually included here, but is essentially just the:

* `get_secondary_structure_DSSP()` function (to calculate helicity)
* `get_radius_of_gyration()` function (to calculate the radius of gyratin), and
* `get_all_SASA()` for calculating solvent-accessible surface area using the `mode='sidechain'` selector.


