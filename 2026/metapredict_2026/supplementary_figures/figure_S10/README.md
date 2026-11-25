# Figure S10 — short IDRs that metapredict V3 misses

Ten CAID3 Disorder-PDB proteins with short IDRs (10–49 residues) that metapredict V3 misses.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A–E | `figure_S10A–E.ipynb` | Loops (IDRs between ordered residues). `figure_S10A.ipynb` also draws the legend and prints the selection. |
| F–H | `figure_S10F–H.ipynb` | IDRs next to unscored residues (neither DisProt-disordered nor resolved in the PDB). |
| I–J | `figure_S10I–J.ipynb` | IDRs at the N- or C-terminus. |
| A–J | `figure_S10.ipynb` | All ten panels with panel letters and one legend (`figures/figure_S10.pdf`): A–E on the left, F–J on the right. |

Each panel has three tracks (the metapredict V3 score with predicted IDRs shaded; AlphaFold2 pLDDT; the Disorder-PDB labels with the missed IDRs outlined) beside the AlphaFold2 model coloured by pLDDT, with the missed IDRs as thick tubes.

`short_idr_examples.py` defines short IDRs, picks the examples and draws the panels. `ground_truth_panels.py` (track drawing, copied from Figure 2), `caid3_panels.py`, `caid_eval.py` and `predictor_classes.py` are the CAID3 code, and `figtools.py` is the shared style. `render_structures.py` downloads the AlphaFold models and renders them with PyMOL.

**Running:** open a notebook from this directory and run all cells (~15 s each). Each writes `figures/figure_S10<panel>.pdf`. Requires metapredict 3.1.x, numpy, matplotlib and scikit-learn. Reads `caid_references/`, `caid3/AlphaFold-pLDDT.caid`, `metapredict_caid/`, `pdb_cache/disprot/` and `alphafold_renders/` from `../../data/` (see its README). The renders are provided; to remake them run `python render_structures.py` (needs `pymol-open-source-whl` 3.2 and, the first time, network access).

**Definitions (Methods)**

- A short IDR is a run of 10–49 residues labelled disordered in CAID3 Disorder-PDB; unscored residues end a run.
- V3 misses an IDR when fewer than half of its residues score at least 0.5 and fewer than half lie in a predicted IDR.
- The examples are the proteins of at most 600 residues with the most missed short-IDR residues. Each context is shown in proportion to its share of all missed short-IDR residues (53% loops, 26% next to unscored residues, 21% tails), and a protein is classed by the context of most of its missed residues.
- The pLDDT track is CAID3's AlphaFold-pLDDT prediction, which is distributed as 1 − pLDDT/100 and converted back.
- Structures are AlphaFold DB models (`AF-<UniProt>-F1`, version 6), whose sequences match the CAID3 sequences, coloured by their own pLDDT. `../../data/alphafold_renders/models.tsv` gives how closely each model's pLDDT matches CAID3's.

**The examples**

| Panel | Protein | Short IDRs V3 misses (residues) |
| --- | --- | --- |
| A | DP04063 (P07275), 575 aa | 23–42, 379–395, 524–546, 566–575 |
| B | DP04132 (Q9NZC2), 230 aa | 131–174, 189–206 |
| C | DP04127 (W1QJK4), 175 aa | 1–19, 47–56, 115–128, 141–158 |
| D | DP04092 (P02710), 461 aa | 205–222, 356–393 |
| E | DP04392 (Q83RJ4), 571 aa | 269–278, 507–540 |
| F | DP03999 (Q8IL97), 200 aa | 126–174 |
| G | DP04185 (C1LGM4), 546 aa | 74–114 |
| H | DP04193 (P12969), 93 aa | 38–74 |
| I | DP04167 (Q8N6T7), 355 aa | 1–25, 62–84 |
| J | DP04128 (P15126), 296 aa | 1–43 |

Across all 319 proteins, 48% of the 162 short IDRs V3 misses have a mean pLDDT above 70, against 12% of the 131 it finds (printed by `figure_S10A.ipynb`).
