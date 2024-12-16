# IDR-omes
##### Last updated 2024-02-04

## Full human IDR-ome 
We often consider the human proteome as the human reference proteome, which typically is around ~24,000 amino acid sequences that correspond to the reviewed protein sequences in UniProt and then corresponding 'canonical' isoform for each of the canonical protein coding genes in the human proteome.

However, many proteins exist as alternative isoforms, and many additional genes are not considered reviewed. 

Using ALBATROSS, we took the COMPLETE set of human proteins (both reviewed and unreviewed, all isoforms), which corresponds to 104,557 different protein-coding sequences and 

1. Predicted all IDRs (142,222 IDRs )
2. Calculated amino acid properties for those IDRs, and
3. Predicted ensemble-average conformational properties for those IDRs

The resulting full IDR-ome (`IDRome_uniprotkb_proteome_UP000005640_2024_02_04.csv`) is provided here as a single large CSV file that can be opened in Excel or easily parsed into your favorite analysis software.

In addition, we include the complete set of protein sequences used here (`uniprotkb_proteome_UP000005640_2024_02_04.fasta.gz`) as a reference.

Note you can re-create this file in a couple of minutes using the [ALBATROSS-colab IDR-ome constructor notebook](https://github.com/holehouse-lab/ALBATROSS-colab)
