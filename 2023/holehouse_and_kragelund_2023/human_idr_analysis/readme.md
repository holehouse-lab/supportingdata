# Generate human IDR-ome annotation
##### Last updated 2023-08-20

### 0. Contact
For any issues with this notebook, protfasta, metapredict, or shephard, please contact [Alex](https://www.holehouselab.com/) via electronic mail.

### 1. Install required programs
The only package needed to generate the data here is the Python-based disorder predictor [metapredict](https://metapredict.readthedocs.io/).

	pip install metapredict

### 2. Aquire the reference human proteome and predict all IDRs.
1. Navigate to the [UniProt human reference proteome](https://www.uniprot.org/proteomes/UP000005640)
2. Download the compressed proteome file (`UP000005640_9606.fasta.gz`)
3. Clean the input file to ensure all sequences consist only of the 20 standard amino acids. This can be done via command-line using `pfasta`, which comes with the `protfasta` package, a library for reading/writing FASTA files and a required dependency for `metapredict`. Running the commmand:


		pfasta --invalid-sequence convert-res -o human_proteome_clean.fasta UP000005640_9606.fasta 
	
	Should reveal the following output

		........................
		pfasta (0.4)
		Please report bugs to:
		https://github.com/holehouse-lab/protfasta
		........................
		
		[INFO]: Reading in the file UP000005640_9606.fasta
		[INFO]: Read in file with 220685 lines
		[INFO]: Parsed file to recover 20586 sequences
		[INFO]: Converted 91 sequences to valid sequences
		[INFO]: Writing new sequence file [human_proteome_clean.fasta]...
	
	
4. Predict all human IDRs with metapredict. 

		metapredict-predict-idrs --mode shephard-domains-uniprot -o shprd_domains_idrs.tsv --verbose human_proteome_clean.fasta	
	Using metapredict V2-FF on an Apple Macbook pro (2022) this takes ~2.5 minutes. If GPU acceleration is available on a Linux machine, this takes ~40 seconds.  
	
	The resulting file (`shprd_domains_idrs.tsv`) is a TSV file that defines each IDR in the human proteome. 
	
	
### 3. Perform analysis 

The proteome can be read and annotated with the IDRs using the Python-based sequence analysis framework [SHEPHARD](https://shephard.readthedocs.io/). SHEPHARD does require some working knowledge of Python, but was designed to be easy to use.


### 4. Make hnRNA1 figure
For the disorder profile figure of hnRNPA1, this can be generated using the `hnRNPA1_figure.ipynb` notebook. 

	
## References
Emenecker, R. J., Griffith, D. & Holehouse, A. S. Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure. Biophys. J. 120, 4312–4319 (2021).

Emenecker, R. J., Griffith, D. & Holehouse, A. S. Metapredict V2: An update to metapredict, a fast, accurate, and easy-to-use predictor of consensus disorder and structure. bioRxiv 2022.06.06.494887 (2022). doi:10.1101/2022.06.06.494887

Ginell, G. M., Flynn, A. J. & Holehouse, A. S. SHEPHARD: a modular and extensible software architecture for analyzing and annotating large protein datasets. Bioinformatics 39, (2023).

## Coda
If this analysis is of interest, pre-computed fully annotated human proteomes are available as Google colab notebooks for easy proteome-wide bioinformatics. 

**Github link**: [https://github.com/holehouse-lab/shephard-colab](https://github.com/holehouse-lab/shephard-colab)

**Colab notebook**: [link](https://colab.research.google.com/drive/1cUHClcA4Fcl-byP_syIDRwycANpO8Na2)
