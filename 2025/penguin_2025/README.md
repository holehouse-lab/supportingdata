# PENGUIN: reproducing the figures from Lotthammer & Holehouse (2025)

#### Last updated 2026-08-28

#### Associated Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22152422.svg)](https://doi.org/10.5281/zenodo.22152422)

This directory contains code and figures for the paper:

> Lotthammer, J. M. & Holehouse, A. S. *Disentangling folding from energetic traps in simulations of disordered proteins.* **J. Chem. Inf. Model.** (2025). doi:[10.1021/acs.jcim.4c02005](https://doi.org/10.1021/acs.jcim.4c02005)

and explains what is being computed, what the panels mean, and what to look for in the result.

## PENGUIN - a method to assess sampling quality

PENGUIN is a computational method for evaluating sampling quality of disordered protein simulations. The method asks two questions of every residue in a set of replicate simulations. 

* **First**, how different is its backbone (φ, ψ) distribution from the same residue in the *excluded-volume (EV) limit* - a reference ensemble in which the chain cannot overlap itself but has no other interactions, and which is therefore the most conformationally heterogeneous ensemble that sequence can have? 
* **Second**, how different is it between independent replicates? 

The first question, on its own, cannot tell genuine structure from an energetic trap: both restrict a residue relative to the EV limit. Adding the second resolves it. Reproducible structure gives a large EV-distance but a small replicate-to-replicate distance; a trap gives a large EV-distance *and* a large replicate-to-replicate distance, because each replicate is stuck somewhere different. Both distances are Hellinger distances between histogrammed dihedral distributions (0 = identical, 1 = disjoint).

For more details, please check out the paper!

## Notebooks

| notebook | paper figure | system | what it shows |
|----------|--------------|--------|---------------|
| [`01_fig1_energy_landscape_schematic.ipynb`](01_fig1_energy_landscape_schematic.ipynb) | Fig 1 | - (synthetic) | cartoon folded vs disordered landscapes |
| [`02_fig2_excluded_volume_landscape.ipynb`](02_fig2_excluded_volume_landscape.ipynb) | Fig 2B | - (synthetic) | the same disordered landscape in the EV limit |
| [`03_fig4_fs_peptide.ipynb`](03_fig4_fs_peptide.ipynb) | Fig 4 | FS-peptide | the baseline PENGUIN summary: transient helix, replicates agree |
| [`04_fig5_fs_peptide_fixed_region.ipynb`](04_fig5_fs_peptide_fixed_region.ipynb) | Fig 5 | FS-peptide, residues 3-10 frozen | a facsimile of a folded domain: EV-distance ≈ 1 but replicates agree perfectly (all-to-all ≈ 0) |
| [`05_fig6_polyGV_PGSK.ipynb`](05_fig6_polyGV_PGSK.ipynb) | Fig 6 | (GV)20-(PGSK)10 | a hydrophobic patch collapsing into replicate-specific traps; normalised distance map; φ/ψ-only variants |
| [`06_fig7_ash1_and_ntd.ipynb`](06_fig7_ash1_and_ntd.ipynb) | Fig 7 | Ash1, SARS-CoV-2 N NTD | biological IDRs sampled efficiently by extensive Monte Carlo |
| [`07_figS22_polyGV_rg_convergence.ipynb`](07_figS22_polyGV_rg_convergence.ipynb) | Fig S22 | (GV)20-(PGSK)10 | global Rg convergence (ANOVA) despite local traps |
| [`08_figS25_gcn4_temperature.ipynb`](08_figS25_gcn4_temperature.ipynb) | Fig S25 | GCN4 fragment, 150-400 K | temperature improves sampling but melts transient helix |
| [`09_figS5-S12_fs_peptide_bin_width.ipynb`](09_figS5-S12_fs_peptide_bin_width.ipynb) | Figs S5-S12 | FS-peptide | sensitivity of the absolute values to the histogram bin width |
| [`10_figS13-S20_synthetic_sequence_bin_width.ipynb`](10_figS13-S20_synthetic_sequence_bin_width.ipynb) | Figs S13-S20 | 40-residue synthetic sequence | the same sweep for a ~30x better-sampled system; uses the precomputed EV reference |
| [`11_figS24_ash1_molecular_dynamics.ipynb`](11_figS24_ash1_molecular_dynamics.ipynb) | Fig S24 | Ash1, long all-atom MD | one trajectory split into segments as pseudo-replicates; uses the precomputed EV reference |

Figure 3 (a flowchart), Figure 2A/C/D, Figure 6C (a molecular rendering) and Figure S24F (SAXS) were made outside Python and have no notebook. Figure S23 (autocorrelation analysis) is not in the paper's public notebooks and is not reproduced.

All figure code lives in [`penguin_helpers.py`](penguin_helpers.py) so each notebook is short: `load_system` discovers replicate trajectories and builds a `SamplingQuality`; `penguin_figure` draws the four-panel summary with the paper's annotations (highlighted regions, dividers, panel letters). The functions are thin wrappers over the public `soursop.sssampling` API, so they double as a template for your own systems. We note that 

## Data

The simulations behind these figures (~5.8 GB) live in [`data/`](data/README.md), one folder per system with `full/` (full-Hamiltonian) and `ev/` (excluded-volume) replicate directories, each holding `__traj.xtc` + `__START.pdb` in the layout `soursop.sstools.find_trajectory_files` expects. It is a curated subset of the paper's simulation archive containing only the trajectories the notebooks use, and it is **excluded from the SOURSOP git repository** because of its size - contact the authors for a copy. `data/README.md` documents every folder (system, replicate counts, frame counts, provenance).

| folder | system | notebooks |
|--------|--------|-----------|
| `data/fs_peptide/` | FS-peptide (FH coil + helical starts, EV) | 03, 04 (EV only), 09 |
| `data/fs_peptide_fixed/` | FS-peptide with residues 3-10 frozen | 04 |
| `data/polyGV_PGSK/` | (GV)20-(PGSK)10 | 05, 07 |
| `data/ash1/` | Ash1 Monte Carlo | 06 |
| `data/ntd/` | SARS-CoV-2 nucleocapsid NTD | 06 |
| `data/gcn4/` | GCN4 fragment at ten temperatures + EV at 270 K | 08 |
| `data/synthetic_sequence/` | the synthetic 40-mer (no EV runs) | 10 |
| `data/ash1_md/` | Ash1 all-atom MD in five segments - **not redistributable** (D. E. Shaw Research; see `data/ash1_md/README.md` for how to obtain it) | 11 |

If you keep the data elsewhere, set `PENGUIN_DATA_DIR` to that directory; `penguin_helpers.get_data_dir()` honours it.

Data are available via the Zenodo repo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22152422.svg)](https://doi.org/10.5281/zenodo.22152422)

## Running

Run the notebooks top to bottom from this directory (they import `penguin_helpers` relative to the working directory). Notebooks 01-07 take seconds to a couple of minutes; 08 and 09 take a minute or two; 10 re-reads ~2 GB per bin width (about six minutes in total) and 11 loads 170,000 all-atom frames. Figures are written under `figures/`.

`original/` holds the two notebooks exactly as published with the paper, and `reference_figures/` the PDFs they produced, for comparison.

## Differences from the published figures

* **Sequence panels.** Panel A of each PENGUIN figure (the sequence cartoon) was drawn in a vector graphics editor; the notebooks print the sequence with a residue ruler instead.
* **Captions vs figures.** The published captions for Figures 5, 6 and 7 describe the Hellinger distances as φ-only, but the code and the published figure files use the joint 2D (φ, ψ) distribution. The notebooks reproduce the 2D figures; notebook 05 also generates the φ- and ψ-only versions.
* **Schematics.** Figures 1 and 2B are random Gaussian landscapes. The originals used the legacy global `numpy` seed; the helper takes an explicit `numpy.random.Generator`, so the bump placement differs while the construction is identical.
* **Precomputed EV reference (notebooks 10 and 11).** Systems without their own EV simulations use SOURSOP's tabulated EV dihedral distributions. SOURSOP 2.0.6 updates the per-residue lookup to address a very minor bug, so absolute values for these two figures differ slightly from the published PDFs. The systems with explicit EV simulations are unaffected.
* **Figures S19/S20.** The published synthetic-sequence bin-width curves clustered tightly and had sharp peaks at residues 17, 25 and 33; both were caused by the old precomputed-EV code (see previous item). The updated curves with the new method are smooth and still less bin-width sensitive than the FS-peptide, although the effect is roughly halved rather than eliminated.
* **Figure S22 statistic.** The published ANOVA p = 0.96 came because we previously analyzd by grouping the $R_g$ values by frame rather than by replicate. Grouped by replicate the test is overwhelmingly significant, which is actually entirely expected given we're looking at 3200 correlated frames per group. Notebook 07 reports both modes of analysis (for completeness), plus the more meaningful comparison of the spread of replicate means (≈ 0.6 Å) with the within-replicate width (≈ 4.6 Å). The figure itself is unchanged.
* **Uncapped chains (Figure S24).** SOURSOP 2.0.6 restricts the dihedral analysis to residues carrying both φ and ψ; on the uncapped Ash1 MD topology that is residues 2-82. Earlier versions paired φ of residue k+1 with ψ of residue k on uncapped chains.
* **Ash1 MD trajectory.** The all-atom trajectory behind Figure S24 is from Robustelli, Piana & Shaw (PNAS 2018) and cannot be redistributed by us; `data/ash1_md/README.md` explains how to request it and rebuild the five-segment layout. Notebook 11 is the only one that needs it.
* **Ash1 MD segments.** The published S24 was actually built from only three of the five trajectory segments. Notebook 11 uses all five.
* **(GV)20-(PGSK)10 distance map.** The caption says n = 10 replicates; the analysis used all 20, as here.
