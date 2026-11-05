# Figure 2B — Predicted disorder across viral lineages

This directory contains the analysis and plotting behind **Figure 2B**, which
shows how the predicted fraction of intrinsic disorder varies across viral
families. Each point is a single viral proteome; families are laid out along the
x-axis and sorted by their mean disorder, so the figure reads left-to-right from
the least- to the most-disordered viral lineages.

The disorder itself is predicted with **metapredict V3**. The heavy lifting —
running metapredict over every viral proteome — happens upstream; what lives here
are the summarised per-proteome disorder values, the taxonomy needed to group
them into families, and the notebook that turns those into the figure.

## Contents

| File | Description |
| --- | --- |
| `Figure2B.ipynb` | The plotting notebook. Loads the per-proteome disorder values, groups proteomes by viral family, filters, and draws the jittered scatter/strip plot. |
| `Fig2B_AFPP_min_10_per_group_no_single_gene_proteomes.pdf` | The rendered figure (AFPP metric, families with ≥10 proteomes, single-gene proteomes excluded). |
| `../../../data/viral_proteomes/viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` | **Input (prediction data).** Per-proteome disorder from metapredict V3. |
| `../../../data/viral_proteomes/viral_proteome_id_to_count.tsv` | **Input.** Number of proteins in each viral proteome. |
| `../../../data/viral_proteomes/viruses_w_lineage.tsv` | **Input.** Full NCBI taxonomic lineage for each viral proteome. |

The three `.tsv` inputs are intermediates of the upstream disorder-prediction pipeline. They live in `data/viral_proteomes/` at the top of the repository (shared with Figure S13; see the README there), and the notebook reads them from there (`DATA_DIR` in its first code cell).

## The two disorder metrics

Every proteome is summarised in two complementary ways, and the notebook can plot
either (set `disorder_vals` in the final cell):

- **AFPP — Average Fraction Per Protein.** The per-protein disorder fraction,
  averaged over all proteins in the proteome. Every protein contributes equally
  regardless of length, so this reflects the disorder content of a typical
  protein. This is the metric used for the committed figure.

- **TFAP — Total Fraction Across Proteome.** The disordered residues summed across
  every protein, divided by the total number of residues in the proteome. This
  weights long proteins more heavily and reflects the bulk disorder content of the
  proteome's coding sequence.

Both are pre-computed and carried in
`viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` as the
`fraction_disorder_per_protein` and `fraction_disorder_total_proteome` columns.

## Input file formats

All three inputs are keyed on the UniProt proteome ID (`UPxxxxxxxxx`):

- `viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` —
  `organism_id`, `fraction_disorder_per_protein` (AFPP),
  `fraction_disorder_total_proteome` (TFAP).
- `viral_proteome_id_to_count.tsv` —
  `Proteome Id`, `Organism Id`, `Protein count`.
- `viruses_w_lineage.tsv` —
  `Proteome Id`, `Organism`, `Organism Id`, `Protein count`,
  `Taxonomic lineage` (comma-separated NCBI lineage, e.g.
  `Viruses, Riboviria, …, Nairoviridae, Orthonairovirus, Dugbe virus`).

## How the notebook works

The notebook runs top to bottom:

1. **Setup and style.** Imports, and matplotlib settings that keep PDF fonts
   editable in Illustrator/Affinity (`pdf.fonttype = 42`). 
   
2. **Jitter helpers.** `generate_uniformly_spaced_points` and
   `jitter_based_on_density` spread the points of each family horizontally in
   proportion to local density, so dense bands of proteomes fan out into a
   readable strip rather than piling onto one line.

3. **Load counts and lineages.** Reads `viral_proteome_id_to_count.tsv` and
   `viruses_w_lineage.tsv` and builds a list of **single-gene** proteomes (proteome
   protein count == 1).

4. **Load disorder, optionally dropping single-gene proteomes.** Reads the disorder
   TSV into `vir_dis_dict`. With `remove_single_gene_proteome = True` (the default,
   and what the committed figure uses) proteomes with a single protein are excluded —
   they add a lot of noise to the plot.

5. **`break_down_into_groups`.** Assigns each proteome to a viral group by matching
   the requested group names against its taxonomic lineage, then drops any group
   with fewer than `min_per_group` proteomes.

6. **`graph_lineages`.** Draws the figure: a density-jittered strip of points per
   group with a short black bar at each group's mean, x-tick labels rotated, and
   optional shaded bands.

7. **Enumerate viral groups.** Pulls the first rank below `Viruses` from every
   lineage (family, or an unclassified/environmental label) as the candidate groups
   to plot.

8. **Assemble and render.** Sets `min_per_group = 10` and `disorder_vals =
   'per_protein'` (AFPP), computes each group's mean, orders groups by mean disorder,
   and saves
   `Fig2B_AFPP_min_10_per_group_no_single_gene_proteomes.pdf`.

## Reproducing the figure

With the three `.tsv` inputs in `data/viral_proteomes/`, open `Figure2B.ipynb` and run all cells.
The metrics and filters are set at the top of the relevant cells:

- `remove_single_gene_proteome` (cell 4) — drop single-protein proteomes.
- `min_per_group` (final cell) — minimum proteomes for a family to appear.
- `disorder_vals` (final cell) — `'per_protein'` (AFPP) or `'total_proteome'`
  (TFAP). Note: this sets only the y-axis label and the `AFPP`/`TFAP` tag baked into
  the output filename. The final cell calls `break_down_into_groups` without `vals`,
  so the plotted values are always AFPP (this is the bug that made the published
  Figure S6B a copy of 2B). The rebuilt Figure S13B (`supplementary_figures/figure_S13`)
  plots TFAP correctly. The final cell also drops the non-taxonomic groups
  (unclassified viruses and environmental samples).

Requires `metapredict`, `numpy`, `matplotlib`, and `scikit-learn`. To regenerate
the `.tsv` inputs from scratch you need to run metapredict V3 over the viral
proteomes in the upstream pipeline; that step is not included here.
