## Additional Supporting Information
###### Last updated 2023-03-08

## Overview
This repository contains additional information for the upcoming manuscript from Boeynaems et al. 2023.

**Aberrant phase separation is a common killing strategy of positively charged peptides in biology and human disease**. 
Boeynaems, S., Rosa Ma, X., Yeong, V., Ginell, G. M., Chen, J.-H., Blum, J. A., Nakayama, L., Sanyal, A., Briner, A., Van Haver, D., Pauwels, J., Ekman, A., Broder Schmidt, H., Sundararajan, K., Porta, L., Lasker, K., Larabell, C., Hayashi, M. A. F., Kundaje, A., Impens, F., Obermeyer, A., Holehouse, A. S. & Gitler, A. D. bioRxiv 2023.03.09.531820 (2023). [doi:10.1101/2023.03.09.531820](https://www.biorxiv.org/content/10.1101/2023.03.09.531820v1)


## Contact
For questions regarding the files in this repository, please contact Alex (alex dot holehouse at wustl dot edu) or Garrett (g dot ginell at wustl dot edu). For any other questions related to the paper, please contact [Steven 'the muscles from Brussels' Boeynaems](https://www.boeynaemslab.org/).

## Contents:

* `/phosphorylation` - Bioinformatics assessing the density of phosphosites in highly abundant positively charged proteins vs. expected distribution by random chance (**Fig. S10** in preprint)


* `/positive_proteins` - all code analyzing abundant, positively charged proteome across various proteomes (**Fig. S9, S11**)

* `/proteome_analysis` - all code for analyzing proteome-wide NCPR/FCR (**Fig. S7**) 

## Software requirements:
#### OS requirements
The bioinformatics analysis in these directories have been tested on 

* **macOS**: Monterey (12.6)
* **Linux**: Ubuntu 18.04, 20.04, 22.04

#### Python dependencies
Our notebooks use standard Python scientific tools:

* numpy - [installation guidelines](https://numpy.org/install/)
* scipy - [installation guidelines](https://scipy.org/install/)
* matplotlib - [installation guidelines](https://matplotlib.org/stable/users/installing/index.html)
* pandas - [installation guidelines](https://pandas.pydata.org/docs/getting_started/install.html)
* seaborn - [installation guidelines](https://seaborn.pydata.org/installing.html)

Additionally, these analyses require several components from the [Holehouse lab](https://www.holehouselab.com/) bioinformatics stack.

#### SHEPHARD
SHEPHARD is a python framework for woring with large protein datasets. [Documentation](https://shephard.readthedocs.io/), [preprint](https://www.biorxiv.org/content/10.1101/2022.09.18.508433v1) and the data required for the analysis here [can be found here](https://github.com/holehouse-lab/supportingdata/tree/master/2022/ginell_2022/shprd_data).

To install SHEPHARD run

	pip install shephard


#### sparrow
sparrow is our in-development sequence analysis package and [can be installed from here](https://github.com/idptools/sparrow). sparrow is currently being developed but the specific features used in these notebooks (NCPR and FCR calculation) are stable and can be used without concern.

