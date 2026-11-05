# metapredict V3 — figure code

Code and data to regenerate every figure in Lotthammer, Hernández-García, Griffith, Weijers, Holehouse & Emenecker, *Metapredict enables accurate disorder prediction across the Tree of Life*. 

## Layout

| Directory | Contents |
| --- | --- |
| `maintext_figures/figure_<N>/` | One directory per main-text figure. |
| `supplementary_figures/figure_S<NN>/` | One directory per supplementary figure. |
| `data/` | Every input the figures and benchmarks read: CAID predictions (downloaded) and reference sets, metapredict 3.1.0 CAID predictions, PDB structural coverage, sequences, AlphaFold models and renders, viral-proteome and tree-of-life tables, and benchmark results. See `data/README.md`. |
| `benchmarks/` | Scripts that time metapredict and write TSVs to `data/benchmarks/`, read by Figures 1E, S5, S6, S7 and S8. Run them on each benchmark machine; see `benchmarks/README.md`. |
| `shared/` | Master copies of the helper modules copied into figure directories (`figtools.py`, `caid3_panels.py`) and `sync_figtools.sh`. |
| `illustrator_files/` | The assembled figures: the Illustrator `.ai` file and its exported PNG for each figure that has one (1, 2, 3, S1, S5, S12, S13, S15, S16). See `illustrator_files/README.md`. |

Each figure directory has a short `README.md`, one notebook per panel (`figure_1A.ipynb`, `figure_S03B.ipynb`, ...) and small helper `.py` files, and writes its panels to `figures/<panel>.pdf`. Figure directories hold code only: every input they read is in `data/` (found through `figtools.shared_data_path`, with a provenance README for each dataset), and the assembled Illustrator figures are in `illustrator_files/`. Notebooks run top to bottom with the figure directory as the working directory.

## Figures

| Figure | Panels |
| --- | --- |
| 1 | A–C CAID1/2/3 MCC rankings; D Msn2 profiles (V1/V2/V3); E human-proteome prediction time |
| 2 | A/B ROC and PR curves vs the CAID3 field; C example traces; D/E best- and worst-predicted proteins with PDB evidence |
| 3 | Disorder across the tree of life (made by co-authors and copied from the original analysis, with data paths updated; Figure 2 in the original submission) |
| S1 | metapredict V1 vs AF2 pLDDT on two proteins |
| S2 | Building the combined disorder/pLDDT training scores |
| S3 | CAID3 rankings by AUC, APS, Fmax and MCC |
| S4 | V1/V2/V3 profiles against the top CAID3 predictors |
| S5 | Throughput, peak RAM, GPU memory and GPU speed-up across versions and hardware |
| S6 | Prediction time vs batch size |
| S7 | Batch vs sequential prediction |
| S8 | `predict_disorder` vs `predict_disorder_stream` memory and throughput |
| S9 | PDB-disorder-monomer benchmark |
| S10 | Short IDRs metapredict V3 misses: ten example proteins, with AlphaFold2 models |
| S11 | CAID3 metrics and rankings by IDR length |
| S12 | Residues correct per 100 on CAID3 |
| S13 | Figure 3 using total fraction across proteome (A tree copied unchanged; B, C) |
| S14 | Screenshot of the metapredict web server (no directory; nothing is generated) |
| S15 | Synthetic and designed proteins |
| S16 | Human transcription factors |

## Requirements

Python ≥ 3.9 with metapredict 3.1.x, numpy, scipy, matplotlib, scikit-learn, protfasta and jupyter; Biopython only to rebuild the PDB cache; psutil (and optionally py-cpuinfo) for the benchmarks; PyMOL (`pymol-open-source-whl`) only to re-render the Figure S10 structures. Figures use Arial and write PDFs with editable text.

## Colours

metapredict V1 `#FF9D33` (orange), V2 `#2978d6` (blue), V3 `#CC2D46` (crimson) in every panel, defined once in `shared/figtools.py`.
