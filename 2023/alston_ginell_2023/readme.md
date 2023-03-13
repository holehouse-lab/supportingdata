## Additional Supporting Information
###### Last updated 2023-02-26

## Overview
This repository contains additional information for the upcoming manuscript from Alston & Ginell et al. 

**Title**: The analytical Flory random coil is a simple-to-use reference model for unfolded and disordered proteins 

**Authors:** Jhullian J. Alston, Garrett M. Ginell, Andrea Soranno, Alex S. Holehouse

**Preprint link**: TBC

## Contact
For questions regarding the files in this repository and the simulations in general, please contact Alex (alex dot holehouse at wustl dot edu).

## Contents:
The content here is pretty self-explanatory. Each figure has its own directory with one or more notebook(s) that generate panels from the paper. In a couple of places, explicit FRC trajectories are also shared.



## Software requirements:
#### OS requirements
The figures were generated on the following systems:

* **macOS**: Monterey (12.6) 
* **Linux**: Ubuntu 22.04

#### Python dependencies
Our notebooks use standard Python scientific tools:

* numpy - [installation guidelines](https://numpy.org/install/)
* scipy - [installation guidelines](https://scipy.org/install/)
* matplotlib - [installation guidelines](https://matplotlib.org/stable/users/installing/index.html)
* pandas - [installation guidelines](https://pandas.pydata.org/docs/getting_started/install.html)

Additionally, you they also require [mdtraj](https://mdtraj.org/), [soursop](https://soursop.readthedocs.io/), and [afrc](https://afrc.readthedocs.io/).

* **mdtraj** can be installed through conda 

		conda install mdtraj # requires conda-forge

* While **soursop** can be installed using `pip` once mdtraj has been installed:

		pip install soursop

* And `afrc` can be installed using `pip` as

		pip install afrc

## References

