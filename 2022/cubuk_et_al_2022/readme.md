## Additional Supporting Information
###### Last updated 2022-08-08

## Overview
This repository contains additional information for the upcoming manuscript from Cubuk et al.

**Title**: TBA

**Authors:** Jasmine Cubuk, Jhullian J. Alston, J. Jeremías Incicco, Alex S. Holehouse, Kathleen B Hall, Melissa Brereton, Andrea Soranno


## Contents:

* `/code` contains a series of subdirectories that contain specific analyses used in the paper.

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

Additionally, you they also require [mdtraj](https://mdtraj.org/), and [soursop](https://soursop.readthedocs.io/) 

* **mdtraj** can be installed through conda 

		conda install mdtraj # requires conda-forge

* While **soursop** can be installed using `pip`:

		pip install git+https://github.com/holehouse-lab/soursop.git


