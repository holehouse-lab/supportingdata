# Figure 2C — Predicted disorder in viruses that infect different host kingdoms

This directory contains the analysis and plotting behind **Figure 2C**, which
compares the predicted fraction of intrinsic disorder in viral proteomes grouped
by the host they infect. Each point is a single viral proteome; the groups are
laid out along the x-axis (prokaryote-infecting, *E. coli*-infecting,
eukaryote-infecting, *H. sapiens*-infecting) with a short black bar at each
group's mean.

The disorder itself is predicted with **metapredict V3**. The heavy lifting —
running metapredict over every viral proteome — happens upstream; what lives here
are the summarised per-proteome disorder values, the host-taxonomy needed to
split proteomes by host kingdom, and the notebook that turns those into the
figure.

## Contents

| File | Description |
| --- | --- |
| `Fig2C.ipynb` | **The figure notebook.** Loads the per-host disorder tables, combines archaea + bacteria into a "prokaryote-infecting" group, draws the density-jittered scatter, and runs the accompanying statistics (Shapiro–Wilk, KS, Mann–Whitney + effect size). |
| `host_to_tax_id.py` | **The upstream classifier.** Reads Virus-Host DB, splits viruses by host kingdom (Eukaryota / Archaea / Bacteria) and by two specific hosts (human = tax 9606, *E. coli* = tax 562), maps each virus to its UniProt proteome, joins in the disorder predictions, and writes the five per-host `*_infecting_viruses.tsv` tables. |
| `fig_infecting_different_orgs.ipynb` | **Exploratory notebook (does not run as-is).** An earlier, wider exploration with extra host groups (*Sulfolobus*, metazoa, vertebrate, mammal, primate). It expects per-host tables that are **not** included in `data/viral_proteomes/`, so it will `FileNotFoundError` on load — kept only to document the exploration. See "Known issues" below. |
| `prokaryote_vs_eukaryote_infecting_viruses_AFPP.pdf` | The rendered figure (AFPP metric). |
| `virushostdb.tsv` | **Input (download), in `data/viral_proteomes/`.** Virus-Host DB dump mapping viruses to their hosts and lineages. |
| `org_id_to_proteome_id.tsv` | **Input (download), in `data/viral_proteomes/`.** UniProt proteome-ID ↔ organism-ID map. |
| `viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` | **Input (prediction data), in `data/viral_proteomes/`.** Per-proteome disorder from metapredict V3. |
| `archaeal_infecting_viruses.tsv`, `bacterial_infecting_viruses.tsv`, `ecoli_infecting_viruses.tsv`, `eukaryotic_infecting_viruses.tsv`, `human_infecting_viruses.tsv` | **Intermediates, in `data/viral_proteomes/`.** Per-host disorder tables written by `host_to_tax_id.py` and read by the notebook. |

Every `.tsv` is either a third-party download or a regenerable intermediate. They live in `data/viral_proteomes/` at the top of the repository (shared with Figure S13; see the README there); the notebooks (`DATA_DIR` in their first code cell) and `host_to_tax_id.py` read and write them there.

## The two disorder metrics

Every proteome is summarised in two complementary ways, and the figure notebook
can plot either (set `disorder_vals` in the plotting cell):

- **AFPP — Average Fraction Per Protein.** The per-protein disorder fraction,
  averaged over all proteins in the proteome. Every protein contributes equally
  regardless of length, so this reflects the disorder content of a typical
  protein. **This is the metric used for the committed figure** and for the
  statistics.

- **TFAP — Total Fraction Across Proteome.** The disordered residues summed across
  every protein, divided by the total number of residues in the proteome. This
  weights long proteins more heavily and reflects the bulk disorder content of the
  proteome's coding sequence.

Both are pre-computed and carried in
`viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` as the
`fraction_disorder_per_protein` (AFPP) and `fraction_disorder_total_proteome`
(TFAP) columns.

## Input file formats

- `virushostdb.tsv` — the Virus-Host DB dump. Tab-separated, 14 columns. The
  columns used here are `virus tax id` (col 0), `host tax id` (col 7), and
  `host lineage` (col 9, a semicolon-separated NCBI lineage whose first token is
  the host kingdom, e.g. `Eukaryota; Opisthokonta; …`).
- `org_id_to_proteome_id.tsv` — `Proteome Id` (`UPxxxxxxxxx`), `Organism Id`.
- `viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` —
  `organism_id` (the UniProt proteome ID), `fraction_disorder_per_protein` (AFPP),
  `fraction_disorder_total_proteome` (TFAP).
- The five `*_infecting_viruses.tsv` tables — `proteome_id`, `per_protein` (AFPP),
  `total_proteome` (TFAP). One row per proteome in that host group.

## How the pipeline works

### 1. Classify viruses by host — `host_to_tax_id.py`

Run this first. It:

1. Reads `virushostdb.tsv` and, for each virus, reads the first token of the host
   lineage to bucket it as **Eukaryota**, **Archaea**, or **Bacteria**, and
   separately flags viruses whose `host tax id` is exactly **9606** (human) or
   **562** (*E. coli*).
2. Reads `org_id_to_proteome_id.tsv` to map each virus's organism ID to its
   UniProt proteome ID.
3. Reads `viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv` to
   attach the AFPP/TFAP disorder values to each proteome.
4. Writes the five per-host tables consumed by the notebook.

Note that the human and *E. coli* groups are **subsets** of the eukaryote and
bacteria groups respectively (a human-infecting virus is also, by lineage, a
eukaryote-infecting virus). This is intentional — they are plotted as their own
finer-grained groups alongside the broad kingdom groups.

### 2. Build the figure — `Fig2C.ipynb`

The notebook runs top to bottom:

1. **Setup and style.** Imports and matplotlib settings that keep PDF fonts
   editable in Illustrator/Affinity (`pdf.fonttype = 42`).
2. **Load the per-host tables** with `read_in_tsv` (keeps only well-formed
   three-column rows).
3. **Jitter helper.** `generate_density_based_scatter_data` uses a Gaussian KDE
   to spread each group's points horizontally in proportion to local density, so
   dense bands fan out into a readable strip.
4. **Combine archaea + bacteria** into a single `prokaryote_infecting` group.
5. **Draw the figure.** For each group it density-jitters the chosen metric,
   scatters the points, and draws a black bar at the group mean, then saves
   `prokaryote_vs_eukaryote_infecting_viruses_<AFPP|TFAP>.pdf`.
6. **Statistics.** Shapiro–Wilk (normality — all groups come out non-normal), then
   the two planned pairwise comparisons (prokaryote vs eukaryote; *E. coli* vs
   *H. sapiens*) with a two-sample **Kolmogorov–Smirnov** test, a
   **Mann–Whitney U** test, and the **common-language effect size** (CLES).

## Reproducing the committed figure

1. The three inputs (`virushostdb.tsv`, `org_id_to_proteome_id.tsv`, `viral_proteome_fraction_disorder_by_prot_and_total_proteome.tsv`) are in `data/viral_proteomes/`.
2. Run `python host_to_tax_id.py` to (re)generate the five per-host tables there.
3. Open `Fig2C.ipynb` and run all cells. It is preset to `disorder_vals =
   'per_protein'` (AFPP), so it writes
   `prokaryote_vs_eukaryote_infecting_viruses_AFPP.pdf` directly; switch to
   `'total_proteome'` if you want the TFAP version instead.

Requires `metapredict`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, and
`seaborn`. Regenerating the disorder `.tsv` from scratch requires running
metapredict V3 over the viral proteomes in the upstream pipeline, which is not
included here.

## Notes and caveats

A few things to be aware of when re-running or extending this analysis (the first
three were latent bugs that have since been fixed in place):

- **AFPP vs TFAP metric (fixed).** The plotting cell used to be left on
  `disorder_vals = 'total_proteome'` (TFAP), which would have saved a
  `..._TFAP.pdf` even though the committed figure and the statistics cells use
  **AFPP** (`per_protein`). It is now set to `'per_protein'` so the notebook
  reproduces `prokaryote_vs_eukaryote_infecting_viruses_AFPP.pdf` out of the box.

- **Jitter seeding (fixed).** The horizontal spread comes from `np.random.normal`
  inside `generate_density_based_scatter_data`, which was unseeded, so the point
  positions drifted on every run. A `np.random.seed(42)` now precedes the plotting
  loop in both notebooks. This is purely cosmetic — the disorder values, the group
  means, and every statistic are unaffected — it just makes the panel reproducible.

- **`host_to_tax_id.py` length guard (fixed).** The row filter was
  `len(i.split('\t')) >= 9`, but the next line read index `[9]` (host lineage),
  which needs at least **10** fields — a re-downloaded dump with a short/trailing-
  empty row would have raised an `IndexError`. The column positions are now
  resolved from the header once, and the guard is `len(fields) > host_lineage_idx`.

- **The human and *E. coli* groups overlap the broad kingdom groups.** By lineage
  a human-infecting virus is also eukaryote-infecting, and an *E. coli*-infecting
  virus is also bacteria-infecting (hence prokaryote-infecting). This is
  intentional — they are finer-grained groups plotted alongside the broad ones —
  but keep it in mind if you sum group sizes.

- **`fig_infecting_different_orgs.ipynb` needs extra data to run fully.** It loads
  five per-host tables (`Sulfolobus_islandicus_infecting_viruses.tsv`,
  `metazoa_infecting_viruses.tsv`, `vertebrate_infecting.tsv`,
  `mammal_infecting.tsv`, `primate_infecting.tsv`) that are **not** included in `data/viral_proteomes/`,
  so the cells that use them (the archaea/bacteria panel and the eukaryote→primate
  panel) will `FileNotFoundError` until those tables are generated. The stale
  `total.set_xticks(...)` calls that previously crashed the single-panel cells have
  been removed, so the final prokaryote-vs-eukaryote cell now runs on the data
  that *is* present. This notebook is kept only as a record of the wider
  exploration; **`Fig2C.ipynb` is the notebook that produces the published panel.**
