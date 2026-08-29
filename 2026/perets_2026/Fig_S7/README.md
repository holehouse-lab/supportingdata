# Figure S7 - predicted hydrodynamic radii for DDX4 NTD variants

STARLING-derived mean $R_h$ for wild-type DDX4 NTD and the three phase-separation-deficient variants, with the DLS value (42 Angstroms) dashed. This is the one-number consequence of Figure S6: every variant that loses long-range intramolecular contacts also expands.

| | $R_h$ (Angstroms) | vs WT |
|---|---|---|
| WT | 41.33 | - |
| 24RtoK | 43.55 | +2.21 |
| 14FtoA | 43.88 | +2.55 |
| 37EDtoA | 43.47 | +2.14 |

## Contents

| path | |
|------|---|
| [`figure_S7_variant_rh.ipynb`](figure_S7_variant_rh.ipynb) | builds the figure |
| `figures/` | output (PDF + PNG) |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib` and [`starling`](https://github.com/idptools/starling); takes seconds.

No local data - it reads four STARLING ensembles (~1.5 GB) from the shared [`../data/`](../data/README.md), which is **not tracked by git**. Set `PERETS_DATA_DIR` if you keep them elsewhere.

## Differences from the published panel

* **The 24RtoK and 14FtoA bars are swapped in the published figure.** The original notebook assembled the values in the order WT, 14FtoA, 24RtoK, 37EDtoA but passed labels in the order WT, 24RtoK, 14FtoA, 37EDtoA, so the two middle bars carry each other's labels. The published panel shows 24RtoK at 43.9 and 14FtoA at 43.6; the correct values are 24RtoK = 43.55 and 14FtoA = 43.88.

  This notebook binds each label to its filename in one table and iterates over it, so they cannot come apart. No conclusion changes - all three variants expand relative to wild-type either way, which is the figure's point.
* The y-axis reads `Rh (A)`; the published panel renders it as $R_h$ in Angstroms, set in the vector graphics editor.
