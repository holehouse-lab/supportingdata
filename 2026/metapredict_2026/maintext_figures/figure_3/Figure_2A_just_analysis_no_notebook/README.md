# Figure 3A — disorder across the tree of life (analysis only)

`GetSppFamilies.R` turns per-proteome disorder into the per-family means and the family-level tree behind Figure 3A (1,480 families from 12,568 cellular proteomes). The tree was then drawn in iTOL and finished by hand, so there is no plotting notebook. This was Figure 2A in the original analysis.

| File | Contents |
| --- | --- |
| `GetSppFamilies.R` | The analysis script (needs `taxize`, `dplyr` and `ape`). |
| `Fig2a.pdf` | The finished panel. |

**What the script does**

1. Reads `UPID_to_TaxID.tsv` (UniProt proteome ID → NCBI taxonomy ID) and fetches each lineage from NCBI with `taxize`, retrying on NCBI's intermittent HTTP 400 errors.
2. Takes the species and family of each proteome and merges in the per-proteome disorder values from `lineages_R.csv`. Viruses and duplicate rows are dropped.
3. Writes the mean and standard deviation per family of both metrics, with the first species as the family's representative: `MeanPerFamily_AFPP.csv` (average fraction per protein, Figure 3A) and `MeanPerFamily_TFAP.csv` (total fraction across the proteome, Figure S13A).
4. Builds a tree from one representative species per family (`class2tree`) and writes it in Newick format. `Families_tree_rearranged.nwk` is that tree after manual rearrangement.

**Data.** The script reads and writes `../../../data/tree_of_life/` (`data_dir` near its top); see the README there. The two CSVs and the rearranged tree are provided. The two inputs (`UPID_to_TaxID.tsv`, `lineages_R.csv`) are large intermediates of the 2024 prediction run and are not distributed.

**Running:** with the two inputs in `data/tree_of_life/`, run `Rscript GetSppFamilies.R` from this directory. It queries the live NCBI Entrez API (set `ENTREZ_KEY` for a higher rate limit), so results can drift as NCBI's taxonomy changes. The files provided are the ones used for the figure.
