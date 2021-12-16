# AlphaFold2 miscellaneous files
##### Last updated 2021-12-16



#### `diff` between V2 and V1 PDB files for human proteome
As a sanity check we ran a `diff` between the PDB files from the first AF2 UniProt dump (V1) and the second AF2 UniProt dump for the human proteome. Specifcally, 


	# run from inside a directory where the v2 files where unpacked
	for d in *pdb 
	do
		d2=$(echo $d | sed "s/v2/v1/g")
		diff ${d} ../../UP000005640_9606_HUMAN/${d2} >> all_diff.txt
	done

The `all_diff.txt` file is here, but spoiler alert is that the following three proteins chained:

* Q04741
* P11586
* Q9Y2G2

And a tiny subset of Arg residues had `NH1` and `NH2` atoms flipped around. 

Other than that no real changes.