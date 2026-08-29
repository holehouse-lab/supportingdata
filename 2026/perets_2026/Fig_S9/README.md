# Figure S9 - inferred hydrodynamic radii for 144 proteins, with DDX4 overlaid

Inferred $R_h$ against chain length for a compilation of disordered proteins with SAXS-measured radii of gyration (red), with the two DDX4<sup>[1-245]</sup> DLS values from this work superimposed at 245 residues: non-phase-separating (blue, 4.2 nm) and phase-separating (yellow, 1.6 nm).

The blue point sits inside the cloud, where a 245-residue disordered protein belongs. The yellow point does not - **only four proteins in the compilation have an inferred $R_h$ below 1.6 nm, and the longest of those is 24 residues.** Of the ten compilation entries within 50 residues of DDX4's length, inferred $R_h$ spans 3.16-4.90 nm. That is the visual statement of the paper's structural claim: on crossing the binodal, DDX4 collapses far below anything a chain of its length ordinarily adopts.

## Contents

| path | |
|------|---|
| [`figure_S9_inferred_rh.ipynb`](figure_S9_inferred_rh.ipynb) | builds the figure |
| [`data/`](data/README.md) | the compiled SAXS radii and sequences, and where they came from |
| `figures/` | output (PDF + PNG) |

## Running

Run the notebook top to bottom from this directory. Needs `numpy` and `matplotlib` only; takes a couple of seconds. Everything it uses is in this directory - it touches none of the bulk simulation data.

$R_g$ is converted to $R_h$ with the Nygaard *et al.* (2017) approximation, with each protein's length taken from its own sequence. That estimator is the consistent choice here because the overlaid DDX4 values are DLS measurements, and Figure S8 shows Nygaard is the estimator that agrees with DLS. As in Figure S2 the conversion is written out rather than imported, since the input is a table of scalar measurements with no coordinates attached.

## Differences from the published panel

* The two DDX4 points are labelled in text with a connecting arrow in the published panel; those were added in a vector graphics editor.
* The original notebook draws the collapsed-state point in yellow, reproduced here. The published panel shows it in dark green, recoloured during figure assembly.
* The published caption says 144 proteins; the CSV holds 145 entries, which is what the notebook plots. Seven pairs share a sequence while being distinct measurements (phosphorylated/unphosphorylated forms, constructs recorded under two names), so the count depends on how those are treated.
* The original computed an AFRC analytical radius of gyration for every sequence and never used it - the slowest step in the original notebook. It has been dropped.
* The original read the CSV with pandas; this uses the standard library `csv` module. Ten entries have a blank error field, which pandas silently made NaN and which is now handled explicitly.
