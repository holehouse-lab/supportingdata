# Scripts
###### Last updated 2021-06-28

## Preamble
This directory contains a simple Jupyter script for computing mean end-to-end distance from the simulation trajectory. 

## Setup
To analyze simulation trajectories you should install the Python-based simulation analysis package [camparitraj](https://camparitraj.readthedocs.io/en/latest/). This can be installed via `pip` and should be relatively painless (note it requires [MDTraj 1.9.5](https://www.mdtraj.org/1.9.5/index.html)).

Alternatively any simulation package will work - we're just computing the end-to-end distance!

## Analysis
See the Jupyter notebook `calculate_end_to_end_distance.ipynb` in this directory