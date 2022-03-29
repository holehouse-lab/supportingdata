## Additional Supporting Information
###### Last updated 2021-12-22

This repository contains additional information for the *in preparation* manuscript:

**The formation of a fuzzy complex in the negative arm regulates the robustness of the circadian clock** <br> Meaghan S. Jankowski, Daniel Griffith, Divya G. Shastry, Jacqueline F. Pelham, Garrett M. Ginell, Joshua Thomas, Pankaj Karande, Alex S. Holehouse, and Jennifer M. Hurley


## Contents:

* Jupyter notebook containing code to recreate certain manuscript figures
* `/data` folder containing various data files used for bioinformatics
* `/simulations` folder containing data and analysis for simulations

For more information on the files in this directory, please contact alex. For information on the remainder of the paper, please contact Jen!

Additional contents can and will be added to this repository as needed. 

## Reproducing analyses & figures:

The sequence and bioinformatic analyses in this manuscript were carried out using Python (version 3.7). To reproduce these analyses and figures download this entire repository, and open the Jupyter Notebook `jankowski_et_al_computational_analyses_and_figures.ipynb`. Then run through the code sequentially, one cell at a time. Depending on how the repository was downloaded, you may need to update some of the filepaths in the notebook to point to the files contained in the `/data` directory.

### Software requirements:
**OS Requirements**
These analyses have been tested on the following systems:
* macOS: Big Sur (11.2.3)
* Linux: Ubuntu 18.04

**Python dependencies**
Analyses use standard Python scientific computing packages:
```
numpy
scipy
pandas
seaborn
matplotlib
```
Additionally, you will also need the *protfasta*, *localCIDER* and *metapredict* packages:
```
pip install protfasta
pip install localcider
pip install metapredict
```
