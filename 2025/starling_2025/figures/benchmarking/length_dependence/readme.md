# Length-dependence benchmarking

We benchmark how STARLING performs as a function of length by sequentially generating ensembles for a 50 amino acid 

The script `benchmark_length_dependence.py` does this benchmarking. This script provides a command-line tool with the following function signature:


	usage: benchmark_length_dependence.py [-h] [--device DEVICE] [--sleepytime SLEEPYTIME]
	
	Run the Starling generate function with adjustable parameters.
	
	options:
	  -h, --help            show this help message and exit
	  --device DEVICE       Device to use for computation (e.g., cpu, cuda, mps)
	  --sleepytime SLEEPYTIME
	                        Time to sleep between runs in seconds
	                        

This can then be run on different hardware to generate the `time_vs_seqlen_<device>.csv` files, which report on sequence length vs. generation time. In general, there is a very weak dependence on sequence length.
	                        