# About
This script makes use two files:

* `4932.protein.aliases.v12.0.txt` - obtained from [StringDB downloads page](https://string-db.org/cgi/download) after filtering for _Saccharomyces Cerevisiae_
* `4932.protein.physical.links.v12.0.txt` - also obtained filtering for _Saccharomyces Cerevisiae_


The script `yeast_interactome.ipynb` takes these files and - along with the yeast proteome FASTA file - builds a SHEPHARD file that relates each UniProt ID in the proteome to 0 or more 'interactor'  UniProt IDs.

Note that the analysis here does NOT, by default, filter by experimental confidence, but this can be done in the `yeast_interactome.ipynb` if needed.