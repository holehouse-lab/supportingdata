# Figure S13A — family-level tree (total fraction across the proteome)

`Supplemental_Fig_6a.pdf` is the finished panel (it keeps its name from the original analysis, where this was Figure S6). It was made like Figure 3A: `GetSppFamilies.R` (taxize) computes the mean TFAP per family and the family tree, which was then drawn in iTOL and finished by hand.

`GetSppFamilies.R` is identical to the copy in `maintext_figures/figure_3/Figure_2A_just_analysis_no_notebook/`; see the README there for what it does and how to run it. Its outputs for this panel, `MeanPerFamily_TFAP.csv` and `Families_tree_rearranged.nwk`, are in `../../../data/tree_of_life/`.
