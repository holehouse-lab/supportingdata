# About

## Zipped and unzipped files
For distribution in our GitHub repository, we provided compressed (zipped) versions of the following files:

* `2025_pLDDT_scores_SHPRD.tsv.zip`
* `2025_disorder_scores_SHPRD.tsv.zip`
* `disorder_scores_SHPRD.tsv.zip`
* `pLDDT_scores_SHPRD.tsv.zip`
* `conservation_scores_SHPRD.tsv.zip`
* `yeast_sequence_dataset.fasta.zip`

To use the code in this repository, these files must first be unzipped. The unzipped version of these files are named in the local `.gitignore` file, which excludes them from inclusion in the git repository.

In our Zenodo repository, we provide the zipped and unzipped versions.


## Versions of data
Between our original analysis and the year of publication, five long years passed. As such, there are a couple of proteome-wide predictions where the precise numerical values generated with metapredict in 2025 (V3) differ slightly from those generated in 2020 during the initial analysis. We note that these differences do not change the overall meaning, but the data and figures in the final manuscript use numerical values calculated with metapredict V1. However, Metapredict V3 can be used to recalculate these values, yielding ever-so-slightly different values. For completeness, we include both sets of values here, although note that the `2025_` were not actually used in this manuscript.


### Disorder scores

* Used in manuscript (calculated using metapredict V1): `disorder_scores_SHPRD.tsv`
* As calculated today with metapredict V3 using V1 backwards compatiblity mode: `2025_disorder_scores_SHPRD.tsv` 


### Predicted pLDDT scores
* Used in manuscript (calculated using metapredict V1): `pLDDT_scores_SHPRD.tsv`
* As calcualted today with metapredict V3 (`2025_pLDDT_scores_SHPRD.tsv`).