## Additional Supporting Information
###### Last updated 2022-03-29

This repository contains additional information for the *under review* manuscript:

**The formation of a fuzzy complex in the negative arm regulates the robustness of the circadian clock.**

Jankowski, M. S., Griffith, D., Shastry, D. G., Pelham, J. F., Ginell, G. M., Thomas, J., Karande, P., Holehouse, A. S., & Hurley, J. M. (2022). In bioRxiv (p. 2022.01.04.474980). [https://doi.org/10.1101/2022.01.04.474980](https://doi.org/10.1101/2022.01.04.474980)


## Contents:

* Jupyter notebook containing code to recreate certain manuscript figures
* `/data` folder containing various data files used for bioinformatics
* `/simulations` folder containing data and analysis for simulations

For more information on the files in this directory, please contact [Alex or Dan](http://holehouse.wustl.edu/). For information on the remainder of the paper, please contact [Jen](https://homepages.rpi.edu/~hurlej2/lab_members.html)!

Additional contents can and will be added to this repository as needed. 

## Reproducing analyses & figures:

The sequence and bioinformatic analyses in this manuscript were carried out using Python (version 3.8). To reproduce these analyses and figures download this entire repository, and open the Jupyter Notebook `jankowski_et_al_computational_analyses_and_figures.ipynb`. Then run through the code sequentially, one cell at a time. Depending on how the repository was downloaded, you may need to update some of the filepaths in the notebook to point to the files contained in the `/data` directory.

### Software requirements:
#### OS requirements
The sequence analyses have been tested on the following systems:

* **macOS**: Big Sur (11.2.3) and Monterey (12.2.1)
* **Linux**: Ubuntu 18.04

#### Python dependencies
Analyses use standard Python scientific computing packages (specific versions are named in the `jankowski_et_al_computational_analyses_and_figures.ipynb` notebook:

* numpy - [installation guidelines](https://numpy.org/install/)
* scipy - [installation guidelines](https://scipy.org/install/)
* pandas - [installation guidelines](https://pandas.pydata.org/docs/getting_started/install.html)
* seaborn - [installation guidelines](https://seaborn.pydata.org/installing.html)
* matplotlib - [installation guidelines](https://matplotlib.org/stable/users/installing/index.html)

Additionally, you will also need the [protfasta](https://protfasta.readthedocs.io/en/latest/), [localCIDER](http://pappulab.github.io/localCIDER/) and [metapredict](https://github.com/idptools/metapredict) (version 2.2) packages. These can be installed using `pip`:

	pip install protfasta
	pip install localcider
	pip install metapredict==2.2

