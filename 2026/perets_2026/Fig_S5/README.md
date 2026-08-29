# Figure S5 - global conformational properties from STARLING ensembles

Two 2D histograms with KDE contours over a 6,000-conformer STARLING ensemble: **(A)** radius of gyration vs asphericity, **(B)** hydrodynamic radius vs asphericity, with the DLS value (42 Angstroms) marked.

This is Figure S4 repeated with a completely different method - a deep-learning sequence-based predictor rather than month-long all-atom Monte Carlo. They agree: STARLING gives mean $R_g$ = 41.7 and mean $R_h$ = 41.4 Angstroms, against 45.9 and 44.0 from the all-atom simulations and 42 from DLS, with the same broad, diagonally elongated distribution.

## Contents

| path | |
|------|---|
| [`figure_S5_starling_global_dimensions.ipynb`](figure_S5_starling_global_dimensions.ipynb) | builds both panels |
| `figures/` | output, one PDF + PNG per panel |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib`, `seaborn`, [`starling`](https://github.com/idptools/starling) and [`soursop`](https://soursop.readthedocs.io/); takes well under a minute.

No local data - it reads the wild-type STARLING ensemble (~1.5 GB) from the shared [`../data/`](../data/README.md), which is **not tracked by git**. Set `PERETS_DATA_DIR` if you keep it elsewhere.

$R_g$ and $R_h$ come from the STARLING ensemble object. Asphericity needs the gyration tensor and so cannot be recovered from a distance map; it is computed with SOURSOP from the reconstructed ensemble trajectory.

## Differences from the published panels

* The published figure was assembled in a vector graphics editor; the A and B letters were added there.
* In the original notebook the asphericity line was commented out and the variable survived in the kernel from an earlier session, so the notebook as saved would not have run from a clean start. Here it is computed explicitly, with an assertion that the conformer counts match.
* The original repeated the plotting code once per panel; here it is one function called twice.
