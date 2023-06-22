## Additional Supporting Information
###### Last updated 2023-06-22

## Overview
This repository contains additional information for the paper from Alston & Ginell et al. 

**Title**: The Analytical Flory Random Coil Is a Simple-to-Use Reference Model for Unfolded and Disordered Proteins

**Authors:** Jhullian J. Alston, Garrett M. Ginell, Andrea Soranno, Alex S. Holehouse

**Preprint link**: [https://www.biorxiv.org/content/10.1101/2023.03.12.531990v2](https://www.biorxiv.org/content/10.1101/2023.03.12.531990v2)

**Paper link**: [https://pubs.acs.org/doi/10.1021/acs.jpcb.3c01619](https://pubs.acs.org/doi/10.1021/acs.jpcb.3c01619)

## Citation
**To cite, please use:**

Alston, J. J., Ginell, G. M., Soranno, A. & Holehouse, A. S. The Analytical Flory Random Coil Is a Simple-to-Use Reference Model for Unfolded and Disordered Proteins. J. Phys. Chem. B 127, 4746–4760 (2023).

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

## Documentation
AFRC documentation is at [https://afrc.readthedocs.io/en/latest/](https://afrc.readthedocs.io/en/latest/)

## References

Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S. SOURSOP: A Python package for the analysis of simulations of intrinsically disordered proteins. bioRxiv (2023). doi:10.1101/2023.02.16.528879



