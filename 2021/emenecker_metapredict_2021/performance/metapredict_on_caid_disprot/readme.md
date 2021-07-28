## Performance for metapredict

To assess how quickly metapredict can compute per-residue disorder scores for the 652 sequences in the CAID DisProt dataset, we ran the following command after installing metapredict:


	time metapredict-predict-disorder ../../data/CAID/CAID_DisProt_Dataset.fasta  
	72.83s user 0.96s system 388% cpu 18.975 total
	
As noted above, this shows that it takes 18.975 seconds to calculate per-residue disorder for each of the sequences in the CAID DisProt dataset.

We also installed a local version of AUCPreD, and compared per-sequence prediction times on our local hardware with the CAID hardware to generate an approximate conversion factor with respect to execution time of ~0.5 (i.e. our hardware is 2x faster than the default CAID hardware).

Based on this qualitative assessment, we estimate metapredict's execution time for the full DisProt dataset on the CAID hardware to be ~40 seconds (2*18.975). 

We emphasize this is a qualitative/ballpark number, and for a _bona fide_ quantitative assessment we wait for the next CAID competition. 

