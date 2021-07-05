# Simulation data

## Files
Simulations were run using the ABSINTH implicit solvent model [1] using the ion parameters of Mao et al. [2].

* `run.key` - keyfile for running all-atom simulations. Note that the data in the paper report ensemble-average results from 20 independent simulations.
* `start.pdb` - starting topology file for simulations
* `traj.xtc` – Subtrajectory for 1000 conformations taken from 20 independent simulations. Provided as representative data. The complete ensemble (90 MBs) is available from Alex on both reasonable and unreasonable request - just trying to avoid filling up our GitHub space quota!
*  `DSSP_H.csv` - ensemble average per-residue helicity (as shown in supplementary figure)
*  `RG.csv` - per-conformation radius of gyration for the full ensemble
*  `Re.csv` - per-conformation end-to-end distance for the full ensemble

## Summary
Total ensemble size:

* 30,000 conformations
* Apparent scaling exponent 0.57


## Refs

[1] Vitalis A, Pappu RV. ABSINTH: A new continuum solvation model for simulations of polypeptides in aqueous solutions. J Comput Chem. 2009;30: 673–699.

[2] Mao AH, Pappu RV. Crystal lattice properties fully determine short-range interaction parameters for alkali and halide ions. J Chem Phys. 2012;137: 064104.


