# Shared helper modules

Master copies of the helpers that are copied into the figure directories, so each directory runs on its own while all of them share one style.

| File | Contents | Copied into |
| --- | --- | --- |
| `figtools.py` | Manuscript style and colours, `save_panel`, `shared_data_path` (finds `data/`) and the metapredict 3.1.x check. | Every figure directory except `figure_3`. |
| `caid3_panels.py` | Loads the CAID3 reference set and predictions, scores predictors, and draws the ranking, curve and trace panels. | `figure_2`, `figure_S03`, `figure_S04`, `figure_S09`, `figure_S10`, `figure_S11`. |
| `sync_figtools.sh` | Copies `figtools.py` into every figure directory (`bash shared/sync_figtools.sh`). | — |

Edit the master copy, then re-copy it. `caid3_panels.py` is copied by hand. `perf_plots.py` (Figures 1E, S5–S8) has no master here: its five copies are kept identical by hand.
