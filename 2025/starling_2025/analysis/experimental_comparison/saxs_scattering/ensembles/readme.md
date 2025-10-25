# Ensembles
This directory holds the .starling ensembles used for SAXS comparison. These ensembles were generated with 600 conformers using the sequences in the `experiment` parent directory (i.e., from here `../experiment/sequences.fasta`). 

Note, we then compute scattering curves using the script:

	build_scattering_profiles.sh
	
This is just a wrapper around `foxs`, and calculates scattering profiles for 400 different conformers and then averages them to give the final scattering curve.	

