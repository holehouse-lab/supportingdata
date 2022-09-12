## Additional Supporting Information
###### Last updated 2022-08-08

## Overview
This repository contains additional information for the upcoming manuscript from Cubuk et al.

**Title**: Disordered N-terminal tail of SARS CoV-2 Nucleocapsid protein forms a dynamic complex with RNA

**Authors:** Jasmine Cubuk, Jhullian J. Alston, J. Jeremías Incicco, Alex S. Holehouse, Kathleen B Hall, Melissa Brereton, Andrea Soranno


## Contents:

* `/code` contains a series of subdirectories that contain specific analyses used in the paper.
Generate_Starting_States - A jupyter notebook with instructions to generate IDR+Folded Domain and RNA starting structures and conformations
RootMeanSquareDistance - A jupyter notebook with instructions on how to calculate root mean square distance when in the bound/unbound state
Example_Binding_Calculator - A jupyter notebook with instructions on how to calculate binding affinities and correct for finite size effects

* `/data` contains either trajectories or other data used by the analysis routines in `/code`.


## Software requirements:
#### OS requirements
The sequence analyses have been tested on the following systems:

* **macOS**: Monterey (12.3.1) and Monterey (12.3.1)
* **Linux**: Ubuntu 18.04

#### Python dependencies
Our notebooks use standard Python scientific tools:

* numpy - [installation guidelines](https://numpy.org/install/)
* scipy - [installation guidelines](https://scipy.org/install/)
* matplotlib - [installation guidelines](https://matplotlib.org/stable/users/installing/index.html)
* pandas - [installation guidelines](https://pandas.pydata.org/docs/getting_started/install.html)
* seaborn - [installation guidelines](https://seaborn.pydata.org/installing.html)


Additionally, you they also require [lammpstools](https://github.com/holehouse-lab/lammpstools), [mdtraj](https://mdtraj.org/), and [soursop](https://soursop.readthedocs.io/) 

* **mdtraj** can be installed through conda 

		conda install mdtraj # requires conda-forge

* While **soursop** can be installed using `pip`:

		pip install git+https://github.com/holehouse-lab/soursop.git


