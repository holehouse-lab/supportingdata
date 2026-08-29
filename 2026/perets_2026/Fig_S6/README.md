# Figure S6 - intramolecular scaling maps for wild-type DDX4 and variants

Four maps of every pairwise residue distance, each normalised by the equivalent distance in a serine-proline dipeptide repeat of matching length - a null model for an expanded, uniformly interacting IDR. Values below the top of the scale mark pairs closer together than chain connectivity alone would put them.

**(A) Wild-type** shows extensive long-range compaction, most prominently between residues 80 and 140, where a high local density of arginine-aromatic contacts sits. **(B) 24RtoK**, **(C) 14FtoA** and **(D) 37EDtoA** each remove one class of interaction (cation-pi; cation-pi and pi-pi; electrostatics) and all three lose those contacts globally - mean normalised distance over well-separated pairs rises from 0.666 in wild-type to 0.72-0.76.

The 24RtoK variant uses the DDX4<sup>[1-236]</sup> construct to match prior experimental work, so it is normalised against a 236-residue SP repeat rather than the 245-residue one.

## Contents

| path | |
|------|---|
| [`figure_S6_scaling_maps.ipynb`](figure_S6_scaling_maps.ipynb) | builds all four panels |
| `figures/` | output, one PDF + PNG per panel |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib` and [`starling`](https://github.com/idptools/starling); takes under a minute.

No local data - it reads six STARLING ensembles (~1.5 GB) from the shared [`../data/`](../data/README.md), which is **not tracked by git**. Set `PERETS_DATA_DIR` if you keep them elsewhere.

The colour scale is fixed at 0.60-0.80 across all four panels; holding it fixed is what makes them comparable.

## Differences from the published panels

* The published figure was assembled in a vector graphics editor; the A-D letters were added there. Panel titles are drawn by the notebook so each panel is self-identifying on its own.
* The original built an `AnalyticalFRC` reference and then never used it - normalisation is against the SP-repeat STARLING ensemble in both versions. The unused object has been dropped.
* The function now raises a clear error on a length mismatch between sequence and null model, which matters for the 236-residue 24RtoK.
* The undefined diagonal (0/0) is explicitly wrapped rather than emitting a bare `RuntimeWarning`; behaviour is unchanged and those cells are drawn white.
* `matplotlib.cm.get_cmap` was removed in Matplotlib 3.11, so this uses `matplotlib.colormaps[...]`. The `plasma` colormap is unchanged.
