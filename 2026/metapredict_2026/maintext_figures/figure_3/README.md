# Figure 3 — disorder across the tree of life

This figure was made by co-authors, and its code is copied here from the original supporting-data directory; the only changes are the bug fixes described in each subdirectory's README and the data paths, which now point at the repository's `data/` directory. It was Figure 2 in the original submission and became Figure 3 when it swapped places with the CAID3 figure (November 2026). The subdirectories and files inside keep the original numbering (`Figure_2A…`, `Figure_2B/`, `Fig2C.ipynb`, ...) so they match the published supporting data. The assembled figure (`figure_3.ai`, `figure_3.png`) is in `../../illustrator_files/figure_3/`.

| Panel | Directory | Contents |
| --- | --- | --- |
| A | `Figure_2A_just_analysis_no_notebook/` | R analysis (taxize) producing the per-family disorder means and the family tree (in `../../data/tree_of_life/`) that was drawn in iTOL; no plotting notebook. |
| B | `Figure_2B/` | Viral proteome disorder by viral realm/family (notebook; inputs in `../../data/viral_proteomes/`). |
| C | `Figure_2C/` | Prokaryote- vs eukaryote-infecting viruses (notebooks and host-splitting script; inputs in `../../data/viral_proteomes/`). |
| — | `tree_assets/` | The tree from panel A without panel lettering (PDF/PNG/SVG), for talks. |

Each subdirectory has its own README. The inputs are shared with Figure S13 (`../../data/viral_proteomes/`, `../../data/tree_of_life/`). The underlying disorder predictions are the 2024 metapredict V3 run over 23,129 UniProt reference proteomes; they were not regenerated with metapredict 3.1.0. The notebooks in `Figure_2B/` and `Figure_2C/` follow the original 2024 style and colours, not the manuscript-wide style used in the other figure directories.
