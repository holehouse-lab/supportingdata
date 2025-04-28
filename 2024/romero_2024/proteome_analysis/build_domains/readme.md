## Build domains
The code in these two notebooks takes either the raw AlphaFold PDB files (`build_dodo_domains.ipynb`) or the output from chainsaw (`build_chainsaw_domains.ipynb`) and generates a series of directories in `../data/domains_<x>` where `<x>` depends on the specific algorithm, and within each directory are one or more PDB files with the format `<UID>_<START>_<END>.pdb`, e.g. `P38149_435_631.pdb`. 

This is the first set of scripts that need to be run for domain decomposition.


### Dependencies

For DODO domain decomposition you need the DODO package, which is available at [https://github.com/idptools/dodo](https://github.com/idptools/dodo).

For CHAINSAW, please use the implementation which is available at [https://github.com/JudeWells/Chainsaw](https://github.com/JudeWells/Chainsaw) and also cite the CHAINSAW paper [1].

### References
[1] Wells, J., Hawkins-Hooker, A., Bordin, N., Sillitoe, I., Paige, B. & Orengo, C. Chainsaw: protein domain segmentation with fully convolutional neural networks. Bioinformatics 40, btae296 (2024).
  


