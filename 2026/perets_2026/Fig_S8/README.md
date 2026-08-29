# Figure S8 - comparing experimental and STARLING-derived hydrodynamic radii

Four bars in two colour-coded pairs, each pairing an experimental $R_h$ measurement with the STARLING prediction computed the corresponding way.

| | experiment | STARLING |
|---|---|---|
| green | PFG-NMR, 31.6 (Brady *et al.* 2017) | Kirkwood-Riseman, **33.10** |
| purple | DLS, 42.0 ± 1.1 (this work) | Nygaard *et al.*, **41.42** |

The two experiments disagree by more than 30%, and the figure explains why that is not a contradiction: **one single ensemble reproduces both**, depending only on which estimator is applied to it - the two differ by 8.3 Angstroms on identical coordinates. Kirkwood-Riseman tracks PFG-NMR, Nygaard tracks DLS. This is why the estimator is stated explicitly wherever an $R_h$ appears in this paper; Figures S2, S4, S5 and S7 all compare against DLS and so all use the Nygaard form.

## Contents

| path | |
|------|---|
| [`figure_S8_rh_comparison.ipynb`](figure_S8_rh_comparison.ipynb) | builds the figure |
| `figures/` | output (PDF + PNG) |

## Running

Run the notebook top to bottom from this directory. Needs `numpy`, `matplotlib` and [`starling`](https://github.com/idptools/starling); takes about ten seconds, nearly all of it the Kirkwood-Riseman sum.

No local data - it reads the wild-type STARLING ensemble (~1.5 GB) from the shared [`../data/`](../data/README.md), which is **not tracked by git**. Set `PERETS_DATA_DIR` if you keep it elsewhere.

## Differences from the published panel

* The published panel was restyled in a vector graphics editor - axis frame, tick styling and label placement differ - but the four bar heights, the colours and the DLS error bar are as produced here.
* The two STARLING values are computed live from the ensemble rather than transcribed.
