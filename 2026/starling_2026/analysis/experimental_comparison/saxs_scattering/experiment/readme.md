# About

This directory contains scattering data and their reanalysis using the (wonderful) MFF of [Riback et al.](https://sosnick.uchicago.edu/SAXSonIDPs/).

This means for every SAXS piece of SAXS data here, we:

1. Verified it was of high quality
2. Re-evaluated it in a CONSISTENT (and in our experience, very robust) way


We feel this sets the highest bar for consistent, high-quality experimental data to compare against. We also hope this is generally useful for the community in the future.

### Files

* `sequences.fasta` - the sequences for all the proteins in this directory (53 proteins)
* `mff_analysis_all.csv` - a CSV file we 53 entries with three columns (name, $R_g$ value, $R_g$ error), as calculated for each scattering data directly using the [Molecular Form Factor](https://sosnick.uchicago.edu/SAXSonIDPs/). The Rg/error values here are summaries; for EACH entry, the full MFF output is in the `mff_out.csv` file in each directory.


### Subdirectories
* Each subdirectory (e.g., `alpha_syn`) contains a set of files:
	* `<filename>.dat` (or something like this), which is the original scattering data. If no `_clean.dat` file (described below) is present, this is the file that will be used for the scattering data.
	* `SASD<ID>.dat` - [optional] - some of the files have the original [SASBDB](https://www.sasbdb.org/) scattering data,
	* `<sequence-name>_clean.dat` - the final, cleaned scattering data in a self-consistent format, which, _if available_ will be the scattering data to use. 
	* `mff_out.csv` - full output from the MFF analysis of the cleaned scattering data.