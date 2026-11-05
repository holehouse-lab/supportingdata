# Figure S10 — Short IDRs that metapredict V3 misses: ten example proteins

| Panel | Notebook | Contents |
| --- | --- | --- |
| A | `figure_S10A.ipynb` | Example 1 (loop). Also prints the selection and, across all CAID3 targets, the AlphaFold2 pLDDT of short IDRs V3 misses vs finds. The legend is drawn here only. |
| B–E | `figure_S10B–E.ipynb` | Examples 2–5 (loops). |
| F–H | `figure_S10F–H.ipynb` | Examples 6–8 (next to unscored residues). |
| I–J | `figure_S10I–J.ipynb` | Examples 9–10 (terminal tails). |
| A–J | `figure_S10.ipynb` | All ten panels in one figure with panel letters and one legend (`figures/figure_S10.pdf`): A–E down the left column, F–J down the right. |

Each panel is one full-length CAID3 Disorder-PDB protein. On the left are three tracks: the metapredict V3 score with metapredict's own IDR calls shaded (0.5 threshold dashed), the AlphaFold2 pLDDT (pLDDT = 70 dashed), and the Disorder-PDB labels (disordered, ordered, unscored), with the short IDRs V3 misses outlined. On the right is the protein's AlphaFold2 model as a cartoon in the standard AlphaFold pLDDT colours (> 90, 70–90, 50–70, < 50), with the short IDRs V3 misses drawn as thick, opaque tubes and the rest of the protein translucent. Each notebook prints every short IDR of its protein with V3's coverage and the mean pLDDT.

`short_idr_examples.py` defines short IDRs, picks the examples and draws the panels. It reuses the track drawing and colours of `ground_truth_panels.py` (copied unchanged from Figure 2, which uses it for 2D/E). `render_structures.py` downloads the AlphaFold models and renders them with PyMOL into `../../data/alphafold_renders/` (the notebooks only read these PNGs). `caid_eval.py` and `predictor_classes.py` are the CAID3 scoring code from the original CAID3 re-analysis (unchanged), `caid3_panels.py` loads the data, and `figtools.py` is the shared style.

**Running:** first `python render_structures.py` from this directory (needs PyMOL, `pymol-open-source-whl` 3.2.0.2, and network access on the first run). It writes 2000 × 2000 pixel ray-traced renders, `<DisProt ID>.png`, and `models.tsv` to `../../data/alphafold_renders/`, and stores the models in `../../data/alphafold_models/`. Then open a notebook and run all cells (~15 s each); each writes `figures/figure_S10<panel>.pdf`, and `figure_S10.ipynb` draws all ten into `figures/figure_S10.pdf`. Reads the CAID3 reference set (`caid_references/`), CAID3's AlphaFold-pLDDT predictions (`caid3/AlphaFold-pLDDT.caid`), the metapredict 3.1.0 V3 predictions (`metapredict_caid/`), the cached DisProt records for each protein's name and UniProt accession (`pdb_cache/disprot/`) and the renders (`alphafold_renders/`) from the shared data directory (`../../data`, see its README). Needs metapredict 3.1.x, numpy, matplotlib and scikit-learn. metapredict warns once that one very short CAID3 target is shorter than its domain-smoothing window; that target is not an example. The PyMOL wheel is built against NumPy 1.x, so `render_structures.py` reads coordinates from the PDB files itself and avoids PyMOL's NumPy-returning calls. The panels embed the renders at full resolution (`imshow(..., interpolation="none")`); without that, matplotlib resamples them to 100 dpi in the PDF.

**Definitions:** a short IDR is a maximal run of 10–49 residues labelled disordered (CAID3 has no IDR shorter than 10). It is a *loop* if ordered residues flank both sides, a *terminal tail* if it reaches either terminus, and otherwise *next to unscored* residues (neither DisProt-disordered nor resolved in the PDB). V3 misses an IDR when fewer than half of its residues score ≥ 0.5 and fewer than half fall in a metapredict IDR call (V3's CAID3 Fmax threshold is 0.51, so this matches the per-residue scoring). The pLDDT track is recovered from CAID3's AlphaFold-pLDDT predictions, which are distributed as 1 − pLDDT/100 (AlphaFold-disorder; Piovesan, Monzon & Tosatto, *Protein Sci* 31:e4466, 2022).

**Selection:** proteins of at most 600 residues with at least one missed short IDR are assigned the context holding most of their missed short-IDR residues and ranked within it by missed short-IDR residues. The top five loop proteins, the top three next-to-unscored proteins and the top two terminal-tail proteins are shown, in proportion to each context's share of the short-IDR residues V3 misses on CAID3 (53%, 26% and 21%). These are the clearest cases of each kind, not random draws.

**Structures:** the canonical AlphaFold DB model (`AF-<UniProt>-F1`, now version 6; the earlier version CAID3 used can no longer be downloaded) is used for each protein. Every model's sequence is identical to the CAID3 sequence (checked). Its pLDDT agrees closely with CAID3's values (median absolute difference 0.02–0.88, Pearson r ≥ 0.991; largest single-residue difference 21.5 in E, and at most 11.4 within a missed short IDR, in B); `data/alphafold_renders/models.tsv` lists these per protein. The structure is coloured by the model's own pLDDT; the track shows CAID3's. Each model is oriented on its core (CA atoms with pLDDT ≥ 50: principal axes on the screen axes), choosing the one of four such orientations that turns the missed short IDRs most toward the viewer; the view is framed on the core and the missed short IDRs, so long low-confidence tails can run out of frame.

**Result with metapredict 3.1.0:**

| Panel | Protein | Short IDRs V3 misses (residues; mean pLDDT) |
| --- | --- | --- |
| A | DP04063 (P07275), Δ1-pyrroline-5-carboxylate dehydrogenase, 575 aa | 23–42 (94), 379–395 (95), 524–546 (97), 566–575 (92) |
| B | DP04132 (Q9NZC2), TREM2, 230 aa | 131–174 (51), 189–206 (72) |
| C | DP04127 (W1QJK4), telomere replication protein EST3, 175 aa | 1–19 (56), 47–56 (79), 115–128 (72), 141–158 (73) |
| D | DP04092 (P02710), acetylcholine receptor subunit α, 461 aa | 205–222 (91), 356–393 (47) |
| E | DP04392 (Q83RJ4), E3 ubiquitin-protein ligase IpaH3, 571 aa | 269–278 (50), 507–540 (87); its terminal tails (1–24, 562–571) do not count as missed |
| F | DP03999 (Q8IL97), bis(5′-adenosyl)-triphosphatase, 200 aa | 126–174 (60) |
| G | DP04185 (C1LGM4), asparagine–tRNA ligase, 546 aa | 74–114 (60) |
| H | DP04193 (P12969), islet amyloid polypeptide, 93 aa | 38–74 (72) |
| I | DP04167 (Q8N6T7), sirtuin-6, 355 aa | 1–25 (84), 62–84 (95) |
| J | DP04128 (P15126), ribosomal protein uL18B, 296 aa | 1–43 (95) |

Across all 319 targets, 48% of the 162 short IDRs V3 misses have a mean pLDDT above 70 (median 67), against 12% of the 131 short IDRs it finds (median 48); this holds for loops, tails and IDRs next to unscored residues alike. So about half of the short IDRs V3 misses are ones AlphaFold2 also confidently predicts as folded (e.g. A, J, and the loop in I), while the rest are low-confidence regions V3 still calls ordered (e.g. the TREM2 stalk in B, the 356–393 loop in D, F and G). The examples over-represent the confident kind: 14 of the 20 missed short IDRs shown have mean pLDDT above 70.

New figure; not part of the original analysis.
