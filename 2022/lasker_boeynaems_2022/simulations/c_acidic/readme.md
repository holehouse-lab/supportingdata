## Summary

**Fig. S7B,C**
Simulations run on c_acid PopZ (full-length).

Sequence below for clarity

	>PopZ C-acidic
	MSDQSQEPTM EEILASIRRI ISEPPIFAVP PPPYLAPPPD PPPPPYPASA
	GFPSVPVTTA PRLPVPEPSI PAEAAEADEP ATLDVEEEED EEDDEDDLED
	EEVAEQLVGV SAASAAASAF GSLSSALLMP KDGRTLEDVV RELLRPLLKE
	WLDQNLPRIV ETKVEEEVQR ISRGRGA


## Contents

`RG.csv` - per-conformation radius of gyration (in Angstroms).

`distance_map.csv` - upper-triangle of average inter-residue CA-CA distances (in Angstroms). Scaling maps (Fig. 7B) are generated as the ratio of the true distance map divided by the AFRC distance map (below).

`AFRC_distance_map.csv` - upper-triangle of average inter-residue CA-CA distances (in Angstroms) expected for a sequence-matched polymer behaving as a *bona fide* Gaussian chain (Analytical Flory Random Coil model, preprint forthcoming...).

`full.pdb` - topology PDB file of full-length PopZ

`subtraj.xtc` - sub-sampled trajectory of 5000 conformations, the full (45,000) conformation ensemble is available on literally any request, reasonable or otherwise [it's just too big for GitHub]. 