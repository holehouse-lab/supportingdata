# Data for Figure S2

Two parallel single-column text files. Line *i* of `rgs.txt` is the mean radius of gyration of the DDX4<sup>[1-245]</sup> ensemble simulated at the temperature on line *i* of `temps.txt`. Both have 78 lines.

| file | contents | units |
|------|----------|-------|
| `rgs.txt` | ensemble-average radius of gyration | Angstroms |
| `temps.txt` | Mpipi-GG model temperature | K |

Three entries in `rgs.txt` are `nan` (at 9 K, 95 K and 145 K); those simulations were not analysed and the notebook drops them.

## Where these numbers came from

The $R_g$ values are the reduced output of coarse-grained molecular dynamics simulations of DDX4<sup>[1-245]</sup>, run with the LAMMPS simulation engine and the Mpipi-GG force field (a derivative of the Mpipi force field). As reported in the paper's *Materials and methods*:

* Initial configurations were generated as random walks, followed by up to 1,000 iterations of energy minimisation, continuing until the maximum force dropped below 1 x 10<sup>-8</sup> kcal/mol/Angstrom.
* Simulations were run in the NVT ensemble at 150 mM NaCl, in cubic boxes with a side length of 500 Angstroms and periodic boundary conditions.
* Temperature was held at a series of different values using a weakly coupled Langevin thermostat, with adjustments applied every 100 ps, and a 20 fs integration timestep.
* Each simulation ran for 12 ms with snapshots every 2 microseconds, giving ensembles of 6,000 conformers per temperature.
* A single simulation was run per temperature. The very long timescale of these simulations, combined with the long gaps between snapshots, means each conformer is effectively independent of the preceding one - which is reflected in how smooth the resulting coil-to-globule curve is.

Radii of gyration were computed from those ensembles using [SOURSOP](https://soursop.readthedocs.io/).

**Temperature is a model parameter, not an experimental temperature.** In the Mpipi-GG model temperature scales the strength of intramolecular interactions relative to thermal energy, so this axis is better read as a solvent-quality axis. Absolute values should not be compared against experimental temperatures.

## What is not included

The LAMMPS input decks, force field files and trajectories are not distributed here - only the per-temperature averages that the figure is built from. Contact the authors if you need the underlying simulation data.

## References

Joseph, J. A., Reinhardt, A., Aguirre, A., Chew, P. Y., Russell, K. O., Espinosa, J. R., Garaizar, A. & Collepardo-Guevara, R. *Physics-driven coarse-grained model for biomolecular phase separation with near-quantitative accuracy.* **Nat. Comput. Sci.** 1, 732-743 (2021). doi:[10.1038/s43588-021-00155-3](https://doi.org/10.1038/s43588-021-00155-3)

Thompson, A. P. *et al.* *LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales.* **Comput. Phys. Commun.** 271, 108171 (2022). doi:[10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171)

Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S. *SOURSOP: a Python package for the analysis of simulations of intrinsically disordered proteins.* **J. Chem. Theory Comput.** 19, 5609-5620 (2023). doi:[10.1021/acs.jctc.3c00190](https://doi.org/10.1021/acs.jctc.3c00190)
