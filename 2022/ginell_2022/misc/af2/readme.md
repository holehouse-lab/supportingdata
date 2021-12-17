# AlphaFold2 miscellaneous files
##### Last updated 2021-12-16



#### `diff` between V2 and V1 PDB files for human proteome
As a sanity check we ran a `diff` between the PDB files from the first AF2 UniProt dump (V1) and the second AF2 UniProt dump for the human proteome. Specifcally, 





	#/bin/zsh
	
	
	# run from inside a directory where the AF2 v2 files were unpacked. The V1 files
	# were unpacked in a directory that is at the relative location of 
	# ../../UP000005640_9606_HUMAN/
	for d in *pdb 
	do
		d2=$(echo $d | sed "s/v2/v1/g")
		diff ${d} ../../UP000005640_9606_HUMAN/${d2} >> all_diff.txt
	done

The `all_diff.txt` file is here, but spoiler alert is that the following three proteins changed:

* Q04741 (Homeobox protein EMX1; EMX1)
* P11586 (C-1-tetrahydrofolate synthase, cytoplasmic, MTHFD1)
* Q9Y2G2 (Caspase recruitment domain-containing protein 8, CARD8)

*Probably* because the actual sequences in UniProt changed.

In addition, a  subset of Arg residues had the order of the hydrogen `NH1` and `NH2` atoms flipped around in the PDB file, which obviously does not have any impact on structure, but just FYI.

Other than that no real changes.