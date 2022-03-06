# Data
###### Last updated 2022-03-06

This directory is where simulation files are located.

## Preamble
The data for this paper are provided in a Zenodo directory: 

[https://zenodo.org/record/6332925](https://zenodo.org/record/6332925)

This is a compressed tarball (`.tgz`) which can be expanded and unpacked into the `/data` directory in this directory structure. After doing this the Jupyter notebooks in the parent directory can be run to reconstruct analysis and figures from the paper.

The data are stored in Zenodo to keep the footprint on GitHub down.


## Sub-sampled trajectories

For the linker trajectories, we provided sub-sampled trajectories via `full.pdb` and `subtraj_500.xtc` files. 

Full trajectories are available on request for the simple reason that they're quite big and loading them all onto Zenodo is really time-consuming. However, as soon as we get a request for the full trajectories we'll upload (and please don't be shy about asking!).

That said, sub-selected trajectories are built by uniformly sampling from the full trajectories and provide a good if more sparse representation of the underlying conformational ensemble.


