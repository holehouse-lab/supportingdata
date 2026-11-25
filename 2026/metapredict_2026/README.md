## Overview (metapredict V3)

##### Last updated 2026-11-24

This repository contains files and data for the preprint by Lotthammer et al.:

Lotthammer, J., Hernández-García, J., Griffith, D., Weijers, D., Holehouse, A. & Emenecker, R. Metapredict enables accurate disorder prediction across the Tree of Life. *bioRxiv* 2024.11. 05.622168 (2026) doi:10.1101/2024.11.05.622168.  (V2)

## Layout

| Directory | In git | Contents |
| --- | --- | --- |
| `maintext_figures/figure_<N>/` | yes | One directory per main-text figure (1–3). |
| `supplementary_figures/figure_S<NN>/` | yes | One directory per supplementary figure (S1–S13, S15, S16; S14 is a web-server screenshot, so no notebook here). |
| `benchmarks/` | yes | Scripts that time metapredict and write the TSVs behind Figures 1E and S5–S8. |
| `shared/` | yes | Master copies of the helper modules that are copied into the figure directories. |
| `data/` | no (Zenodo) | Every input the figures and benchmarks read. See `data/README.md`. |
| `illustrator_files/` | no (Zenodo) | The assembled figures (Illustrator files and exported images). |

`data/` and `illustrator_files/` are too large for git and are in the Zenodo archive that accompanies the paper. Unpack `data/` at the top of this directory, or point `$METAPREDICT_V3_DATA_DIR` at it. A notebook that cannot find an input stops with a `MISSING:` message naming the file.

Each figure directory has a `README.md`, one notebook per panel (`figure_1A.ipynb`, `figure_S03B.ipynb`, ...) and small helper `.py` files, and writes each panel to `figures/<panel>.pdf`. Figure directories hold code only. Run notebooks top to bottom with the figure directory as the working directory. Panel letters, figure numbers, and the numbers quoted in each README follow the manuscript.

## Figures

| Figure | Panels |
| --- | --- |
| 1 | A–C CAID1/2/3 MCC rankings; D Msn2 profiles (V1/V2/V3); E human-proteome prediction time |
| 2 | A/B ROC and precision–recall curves against the CAID3 field; C example traces; D best- (left) and worst-predicted (right) proteins with PDB evidence |
| 3 | Disorder across the tree of life (co-author code, copied from the original analysis) |
| S1 | metapredict V1 vs AF2 pLDDT on two proteins |
| S2 | Building the combined disorder/pLDDT training scores |
| S3 | CAID3 rankings by AUC, APS, Fmax and MCC |
| S4 | V1/V2/V3 profiles against the top CAID3 predictors |
| S5 | Throughput, peak RAM, GPU memory and GPU speed-up across versions and hardware |
| S6 | Prediction time vs batch size |
| S7 | Batch vs sequential prediction |
| S8 | `predict_disorder` vs `predict_disorder_stream` memory and throughput |
| S9 | PDB-disorder-monomer benchmark |
| S10 | Short IDRs metapredict V3 misses: ten example proteins with AlphaFold2 models |
| S11 | CAID3 metrics and rankings by IDR length |
| S12 | Residues correct per 100 on CAID3 |
| S13 | Figure 3 using the total fraction across the proteome |
| S14 | Screenshot of the metapredict web server (no directory; nothing is generated) |
| S15 | Synthetic and designed proteins |
| S16 | Human transcription factors |

## Requirements

Python ≥ 3.9 with metapredict 3.1.x, numpy, scipy, matplotlib, scikit-learn, protfasta, and Jupyter. Optional: Biopython to rebuild the PDB cache; psutil (and py-cpuinfo) for the benchmarks; PyMOL (`pymol-open-source-whl`) to re-render the Figure S10 structures; R with taxize, dplyr and ape for the Figure 3A analysis. Figures use Arial and are written in PDFs with editable text.

Figures 1, 2, S1–S12, S15 and S16 use metapredict 3.1.0. Figures 3 and S13 use the 2024 metapredict V3 predictions of 23,129 UniProt reference proteomes, which are not regenerated here.

## Colours

metapredict V1 `#FF9D33` (orange), V2 `#2978d6` (blue), V3 `#CC2D46` (crimson) in every panel, defined once in `shared/figtools.py`.
