# Figure S2 - coil-to-globule curve for DDX4<sup>[1-245]</sup>

Hydrodynamic radius vs temperature from coarse-grained (LAMMPS, Mpipi-GG) simulations, with the 298 K value marked in red. The chain expands ~2-fold, from $R_h \approx 2.2$ nm at low temperature to $\approx 4.5$ nm at high temperature.

Temperature here is a model parameter that scales intramolecular interaction strength, not an experimental temperature - read the x-axis as solvent quality.

## Contents

| path | |
|------|---|
| [`figure_S2_coil_to_globule.ipynb`](figure_S2_coil_to_globule.ipynb) | builds the figure |
| [`data/`](data/README.md) | per-temperature mean $R_g$ values, and where they came from |
| `figures/` | output (PDF + PNG) |

## Running

Run the notebook top to bottom from this directory. Needs `numpy` and `matplotlib` only; takes a couple of seconds. It loads 78 (temperature, $\langle R_g \rangle$) pairs, drops three with no $R_g$ value, converts to $R_h$ using the Nygaard *et al.* (2017) approximation, and plots.

The conversion is written out in the notebook rather than imported. SOURSOP implements the identical arithmetic as `SSProtein.get_hydrodynamic_radius(mode='nygaard')`, but that is a method on a trajectory object, and the input here is a file of scalar averages.

## Differences from the published panel

* The published version was restyled in a vector graphics editor (larger labels, heavier axes, larger marker); the curve and data are identical; this notebook reproduces the raw export.
* The "$R_H$ at 298 K" label was added by hand to the published figure; the notebook draws it with `ax.annotate`.
* Arial now heads a font fallback list so the notebook runs without warnings where Arial is missing.
* The original sliced `[1:]` after dropping NaNs, discarding the 14 K point ($R_h$ = 2.26 nm) as well, so the panel starts at 50 K. Preserved here to match.
