# Figure S4 - global conformational properties of DDX4<sup>[1-245]</sup>

Two 2D histograms with KDE contours describing how chain size and shape co-vary across the all-atom CAMPARI/ABSINTH ensemble: **(A)** radius of gyration vs asphericity, **(B)** hydrodynamic radius vs asphericity, with the DLS measurement (42 Angstroms) marked.

Both show a single broad, diagonally elongated cloud - larger conformations are more elongated - which is what a well-sampled disordered chain looks like. Mean $R_g$ = 45.9 Angstroms and mean $R_h$ = 44.0 Angstroms, the latter in reasonable agreement with the experimental 42 Angstroms, which is what licenses using these simulations to interpret the DLS data. Figure S5 is the same analysis on STARLING ensembles.

## Contents

| path | |
|------|---|
| [`figure_S4_global_dimensions.ipynb`](figure_S4_global_dimensions.ipynb) | builds both panels |
| `figures/` | output, one PDF + PNG per panel |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib`, `seaborn` and [`soursop`](https://soursop.readthedocs.io/) (2.0.5+); takes well under a minute.

This figure has no data files of its own - it reads the CAMPARI DDX4 NTD trajectory (~350 MB) from the shared `data/` directory at the repository root, which is **not tracked by git**. See [`../data/README.md`](../data/README.md); set `PERETS_DATA_DIR` if you keep it elsewhere.

$R_h$ uses SOURSOP's default Nygaard *et al.* estimator - the same approximation as Figure S2, and the one Figure S8 shows agrees with DLS.

## Differences from the published panels

* The published figure was assembled in a vector graphics editor; the A and B letters were added there.
* The published panel B shows the DLS value as a dashed line labelled above the axes. The raw export, and this notebook, draw a solid line with an in-axes legend.
* The trajectory here holds 25,000 conformers, matching this figure's published caption. Note the *Materials and methods* text quotes a final ensemble of 23,773 conformers.
* The original repeated the plotting code once per panel; here it is one function called twice, so the panels are guaranteed to match.
