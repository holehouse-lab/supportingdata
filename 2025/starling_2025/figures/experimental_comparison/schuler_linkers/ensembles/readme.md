# Ensembles directory
On GitHub, `.xtc`, `.pdb`, and `.starling` files are added to `.gitignore` and excluded from the repository; however, these files are included in the Zenodo full repository.


## Ensemble generation
Ensemble here were generated using the following command:

	starling ../schuler_seqs.fasta -c 2000 -r
	
	
## SAXS back calculation
We also back-calculated scattering curves for 

	dCh_minus/
	dTRBP/
	sNh_plus/
	sNrich/
	
Using FOXS. The resulting average scattering curves in each of these subdirectories as `average_curve.dat`. There are additional files in here:

* `average_curve_error.dat` - includes an analytical estimation of the q-dependent error
* `average_curve_noise.dat` - includes stochastic noise from the analytical q-dependent error, and the q-dependent error itself
* `logfile_autoSCRT.txt` - logfile 
		