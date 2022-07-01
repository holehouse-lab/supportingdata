## Additional Supporting Information
###### Last updated 2022-07-01

## Overview
This repository contains additional information for the manuscript:

**The kinetic landscape of human transcription factors**

Nicholas E Mamrak<sup>1</sup>, Nader Alerasool<sup>2</sup>, Daniel Griffith<sup>3</sup>, Alex S Holehouse<sup>3</sup>, Mikko Taipale<sup>2</sup>, Timothée Lionnet<sup>1</sup>*


1. Institute for Systems Genetics, New York University School of Medicine, New York, NY, United States; Department of Cell Biology, New York University School of Medicine, New York, NY, United States

2. Department of Molecular Genetics, University of Toronto, Toronto, ON, Canada; Donnelly Center for Cellular and Biomolecular Research, University of Toronto, Toronto, ON, Canada

3. Department of Biochemistry and Molecular Biophysics, Washington University School of Medicine, St. Louis, MO, United States

4. Center for Science and Engineering Living Systems (CSELS), Washington University in St. Louis, St. Louis, MO, United States


## Contents:

* `/code` contains sequence, preparatory, and analysis code for the transcription factor bioinformatic analysis performed in this paper, and some bonus stuff.

For more information on the files in this directory, please contact [Alex or Dan](http://holehouse.wustl.edu/). For information on the remainder of the paper, please contact [Tim](http://www.timotheelionnet.net/)!

Additional contents can and will be added to this repository as needed. 

## Reproducing analyses & figures:

Parts of the bioinformatic analyses in this manuscript were carried out using the notebooks in `/code` with Python (version 3.7 or 3.8). To reproduce these analyses and figures download this entire repository, and open the Jupyter Notebook `compute_idr_params.ipynb`. 

The code should run directly assuming the required packages are included. For a more detailed discussion of setting up a working conda environment see `code/readme.md`.


### Software requirements:
#### OS requirements
The sequence analyses have been tested on the following systems:

* **macOS**: Big Sur (11.2.3) and Monterey (12.3.1)
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

