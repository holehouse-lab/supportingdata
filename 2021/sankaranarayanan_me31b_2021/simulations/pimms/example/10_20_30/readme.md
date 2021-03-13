# Readme
##### Last updated March 13th 2021

The files in this directory provide ready-to-run example keyfiles for launching simulations with or without RNA at with protein-protein interaction strengths tuned to the weakest value (as shown in Fig. XX).

Specifically these can be run with PIMMS using


	pimms -k KEYFILE_with_RNA.kf
	


#### Interaction strength

The `10_20_30` nomenclature used reflects the strength of IDR:IDR, IDR:OD, and OD:OD interactions, where IDR = N or C terminal disordered beads (beads 1,2,4 and 5) and OD is the ordered domain bead (bead 3). The chain topology is defined in the keyfile as 

	CHAIN : 800 IIOII
	
i.e. these simulations will contain 800 'protein' molecules with the structure `IIOII`.
