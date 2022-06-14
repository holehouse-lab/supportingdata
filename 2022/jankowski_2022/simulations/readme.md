## Simulation information for Jankowski et al.
###### Last updated 2022-03-29

The two subdirectories here contain information from the simulations (`/data`) and to enable re-running the simulations (`/input_files`).


## Software
Simulations were run using version 2.0 of the CAMPARI simulation engine, although equivalent results are obtained running using version 3.0 and version 4.0.

CAMPARI can be installed by following the installation instructions provided  here: [http://campari.sourceforge.net/](http://campari.sourceforge.net/).

Once installed, simulations can be run using the keyfiles in the subdirectories simply using

	campari -k keyfile.kf
	
Assuming the `campari` binary is in the user's PATH. 

Expected run time for simulations is on the order 3-4 days using standard hardware (single CPU, ~1 GB RAM).

The expected output is a collection of files which includes an .xtc trajectory file.
