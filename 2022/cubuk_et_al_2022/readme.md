## Additional Supporting Information
###### Last updated 2022-08-08

## Overview
This repository contains additional information for the upcoming manuscript from Cubuk et al.

**Title**: TBC

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
Analyses use standard Python scientific computing packages (specific versions are named in the `compute_idr_params.ipynb` notebook:

* numpy - [installation guidelines](https://numpy.org/install/)
* scipy - [installation guidelines](https://scipy.org/install/)
* matplotlib - [installation guidelines](https://matplotlib.org/stable/users/installing/index.html)

Additionally, you will also need the cython, [sparrow](https://github.com/holehouse-lab/sparrow), [metapredict](https://github.com/idptools/metapredict) and [shephard](https://shephard.readthedocs.io/) packages. These packages can be installed using `pip`:

	pip install cython
	pip install metapredict
	pip install shephard
	pip install git+https://github.com/holehouse-lab/sparrow.git

As mentioned, for a more detailed discussion of setting up a working conda environment see `code/readme.md`.

