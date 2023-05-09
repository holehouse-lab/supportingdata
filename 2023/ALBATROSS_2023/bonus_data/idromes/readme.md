# IDR-omes
##### Last updated 2023-05-08

The files here are FASTA files with that contain all predicted IDRS in a proteome for a given organism.

In general, we recommend interacting with large-scale IDR annotations as [Domains in SHEPHARD files](http://dx.doi.org/10.1101/2022.09.18.508433).

However, the [ABATROSS colab notebook](https://colab.research.google.com/github/holehouse-lab/ALBATROSS-colab/blob/main/example_notebooks/polymer_property_predictors.ipynb) provides some nice tools for interactively exploring large datasets, so we provide the IDRs here so you can predict sequence-derived properties using these files and then explore in your browser.

If there's a specific organism or set of IDRs you like to see, please let us know.

### Contents
#### 1. human_idrome.fasta
Human IDR-ome with 33,563 IDRs (as predicted using [metapredict V2](http://dx.doi.org/10.1101/2022.06.06.494887). In the colab predicting all radii of gyration takes ~40 seconds.