# Data
###### Last updated 2022-05-25

The files in this directory contain the raw data used for analysis.

Specifically:

* `real_data.tsv` contains the data sent over on activity and intensity. Note this also has a couple of gene names updated (see readme in parent directory):

        CCDC101 ==> CCDC101/SGF29
        KDELC2 ==> KDELC2/POGLUT3
        SDCCAG3 ==> SDCCAG3/ENTR1

* `uniprot_sequences.fasta` is the protein sequences which were obtained by using the set of gene names with UniProt's [ID mapping tool](https://www.uniprot.org/uploadlists/). 

* Note that VP64 is not a real protein, but a 4x tandem of the V16 activation domain (`DALDDFDLDML`). The sequence used is based on info [on the igem website](http://parts.igem.org/Part:BBa_J176013). Note that, slightly oddly, we *don't actually predict VP64 to be disordered...*. Probably not right, but, it's also so short I don't think as a datapoint in this dataset it's necessarily easily comparable with other proteins.