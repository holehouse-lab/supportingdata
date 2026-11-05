# Figure 2A — Predicted disorder across the tree of life

This directory contains the analysis behind **Figure 2A**, which shows how the
predicted fraction of intrinsic disorder varies across families spanning the
tree of life. The workflow takes a set of proteomes, resolves each one to its
NCBI taxonomy, summarises the disorder content per taxonomic family, and builds
a family-level phylogenetic tree that the final figure is drawn onto.

Unlike the other figure directories in this repository, this one ships only the analysis code — there is no plotting notebook here (hence the `_just_analysis_no_notebook` suffix). The tree and the per-family summaries it produced are in `data/tree_of_life/` at the top of the repository; they are the inputs used to render the published figure.

## Contents

| File | Description |
| --- | --- |
| `GetSppFamilies.R` | The analysis script. Resolves taxonomy from NCBI, aggregates disorder per family, and builds the family tree. |
| `../../../data/tree_of_life/MeanPerFamily_TFAP.csv` | Mean disorder per family using the **total fraction across proteome** metric. |
| `../../../data/tree_of_life/MeanPerFamily_AFPP.csv` | Mean disorder per family using the **average fraction per protein** metric. |
| `../../../data/tree_of_life/Families_tree_rearranged.nwk` | The final family-level phylogeny (Newick), one representative tip per family, after manual rearrangement. |
| `Fig2a.pdf` | The rendered figure. |

## The two disorder metrics

Each proteome is summarised in two complementary ways, and both are carried
through to their own per-family CSV:

- **TFAP — Total Fraction Across Proteome.** The disordered residues summed
  across every protein in the proteome, divided by the total number of residues
  in the proteome. This weights large proteins more heavily and reflects the
  bulk disorder content of the organism's protein-coding sequence.

- **AFPP — Average Fraction Per Protein.** The per-protein disorder fraction,
  averaged over all proteins. Every protein contributes equally regardless of
  length, so this reflects the disorder content of a typical protein.

Reporting both separates "how much of the proteome is disordered" from "how
disordered is an average protein", which can diverge when disorder content
correlates with protein length.

## How `GetSppFamilies.R` works

The script runs in five stages. It depends on `taxize` (for NCBI taxonomy),
`dplyr` (aggregation), and `ape` (tree writing).

### 1. Resolve taxonomy from NCBI

The script reads `UPID_to_TaxID.tsv`, a two-column table mapping each UniProt
proteome ID (`V1`) to its NCBI taxonomy ID (`V2`). For every taxonomy ID it
calls `classification(..., db = "ncbi")` to pull the full lineage.

NCBI's Entrez interface intermittently returns `Bad Request (HTTP 400)` errors
that are unrelated to the specific query — a known issue with these bulk
lookups. To work around this the retrieval is wrapped in a loop that catches the
error, waits ten seconds, and retries, so a transient failure does not abort the
whole run.

### 2. Extract species and family per proteome

For each resolved lineage the script pulls out the rows at rank `species` and
rank `family` and pairs them. Lineages that failed to resolve (stored in the
list with length 1) are skipped, and lineages with no family assignment are
padded with an `"empty"` placeholder so the row bind does not fail. The result
is a table of species ↔ family assignments.

### 3. Attach the disorder measurements

The species/family table is merged back onto:

- `UPID_to_TaxID.tsv` — to recover the UniProt proteome ID, and
- `lineages_R.csv` — the per-proteome disorder measurements
  (`Total_fraction_across_proteome` and `Average_fraction_per_protein`).

Viruses are then dropped (any row containing the string `"virus"` is removed)
and duplicate rows are collapsed.

### 4. Summarise per family

The merged table is grouped by family. For each family the script computes the
mean (and standard deviation) of each disorder metric and keeps the first
species as the family's representative. These summaries are written to `data/tree_of_life/` as:

- `MeanPerFamily_TFAP.csv`
- `MeanPerFamily_AFPP.csv`

Each file has one row per family (1480 families here) with columns
`family`, `mean_<metric>`, `sd_<metric>`, and `spp`. The standard deviation is
`NA` for families represented by a single proteome.

### 5. Build the family tree

Using the representative species for each family, the script queries NCBI once
more, converts the classifications to a tree with `class2tree()`, extracts the
`phylo` object, and writes it to Newick with `ape::write.tree()`. This tree —
after manual rearrangement — is provided as `data/tree_of_life/Families_tree_rearranged.nwk`
and defines the branching layout of Figure 2A.

## Reproducing the analysis

The script reads its inputs from, and writes its outputs to, `data/tree_of_life/` at the top of the repository (`data_dir` near the top of the script; the species-tree PDFs it draws are written to this directory). Its two input files are **not** included (they are large intermediate products of the upstream disorder-prediction pipeline):

- `UPID_to_TaxID.tsv` — UniProt proteome ID → NCBI taxonomy ID.
- `lineages_R.csv` — per-proteome disorder metrics, keyed by proteome ID.

With those in `data/tree_of_life/`, run from this directory:

```r
Rscript GetSppFamilies.R
```

Note that the script calls out to the live NCBI Entrez API, so a run requires
network access and its runtime depends on NCBI's responsiveness. Setting an
`ENTREZ_KEY` (via `Sys.setenv(ENTREZ_KEY = ...)`) raises the request rate limit
and makes the bulk lookups more reliable. Because it hits a live service, exact
outputs may vary slightly over time as NCBI's taxonomy is updated — the CSVs and
tree in `data/tree_of_life/` are the versions used for the published figure.  
