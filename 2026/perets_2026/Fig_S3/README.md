# Figure S3 - local conformational properties of DDX4<sup>[1-245]</sup>

Three per-residue panels from all-atom CAMPARI/ABSINTH simulations, asking whether the chain behaves like a random coil locally. **(A)** DSSP helical (blue, up) and extended (red, down) occupancy, violins across 20 replicas. **(B)** Hellinger distance between each residue's phi (purple) / psi (green) dihedral distribution and the same residue in a self-avoiding random coil (SARC). **(C)** Radius of gyration of a sliding 10-residue window, divided by the SARC equivalent.

DDX4<sup>[1-245]</sup> is disordered - no persistent secondary structure anywhere - but is not a featureless random coil: several regions carry transient helicity, dihedrals deviate from SARC along the whole sequence, and local segments are consistently more compact (mean ratio 0.90) than the reference. All three panels agree on the strongest feature, around residues 171-183.

## Contents

| path | |
|------|---|
| [`figure_S3_local_conformational_properties.ipynb`](figure_S3_local_conformational_properties.ipynb) | builds all three panels |
| [`data/`](data/README.md) | per-replica DSSP occupancies (panel A), and provenance |
| `figures/` | output, one PDF + PNG per panel |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib`, `tqdm` and [`soursop`](https://soursop.readthedocs.io/) (2.0.5+); takes about a minute, most of it loading trajectories.

Panels B and C need the two CAMPARI trajectories (~750 MB), which live in the shared `data/` directory at the repository root and are **not tracked by git** - see [`../data/README.md`](../data/README.md). Set `PERETS_DATA_DIR` if you keep them elsewhere. Panel A needs only the DSSP files in `data/`, which are included.

## Differences from the published panels

* The published figure was assembled in a vector graphics editor: panel letters, the shared residue axis, y-axis labels, the "Helicity"/"beta-sheet" and "phi/psi angle" text, and the "Less RC"/"More RC" arrows were all added there.
* The original notebook had a `SHOW_SUBREGIONS` flag shading four helical regions (11-18, 171-183, 195-203, 215-222) in pale blue. The published panel is the unshaded version, reproduced here.
* Panel C window count is taken from the trajectory rather than a hard-coded copy of the sequence; equivalent (236 windows), but removes the duplicated sequence string.
* The original also computed a split-half "self" Hellinger distance as a sampling control (`self_hellinger.pdf`). It is not in the published figure and is not reproduced.
