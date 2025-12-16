# Langstein-Skora, Schmid, and Huth et al.

## Overview
This GitHub repository contains all the data and most of the code associated with generating the figures and performing the analyses associated with the manuscript:

***Sequence and chemical specificity define the functional landscape of intrinsically disordered regions***

written by Langstein-Skora, Schmid, and Huth et al.

The full citation will be updated here when the paper is formally accepted. 

## Changelog
Please see the `changelog_github.md` and `changelog_zenodo.md` for a record of changes made to the two. Changes to Zenodo will lag GitHub.

## Repository structure
For each figure, we provide links to either the Jupyter notebook used to generate the figure from scratch or, for a few cases, the raw data associated with the figures. Given the extent of the data associated with this paper, this directory maps each figure panel to a specific location where data and/or code are provided. However, within each directory, we often also provide README files to help orient the reader and provide context for the files.

In addition to the code in this GitHub repository, ODM-seq analysis was performed using the code in [https://github.com/gerland-group/LS_ODM_seq_analysis](https://github.com/gerland-group/LS_ODM_seq_analysis). For completeness, we provide a mirror of this code in `odmseq_analysis/` directory.

## Source data for publication
In compliance with the journal's requirements, this repository contains all the source data associated with each main text and extended data figure, provided in a way we _hope_ is maximally useful to readers. However, we note that raw numerical values for some figures may not be explicitly saved in this repository if the code used to generate the figures is provided. In principle, to access those numerical values, it is sufficient to re-run the code, and we have checked that all the code here runs in the associated notebooks. *However*, if you run into issues with this and/or running the code is not possible for whatever reason, please contact Alex, and he can provide you with whatever data you might want.

Similarly, if anything about the actual analysis is confusing or you'd like specific data shared in a more convenient format, please contact Alex directly, and he can update the code accordingly to provide this!

## Zenodo vs. GitHub
To enable easy responsiveness and code access, we have made this repository available on GitHub. However, the full size of this repo is ~18 GB (!), which is obviously too big for GitHub. As such, we use `.gitignore` files to exclude large raw data files or output directories with large amounts of generated data. The FULL repository (including all data) is provided at: [10.5281/zenodo.17770967](https://zenodo.org/records/17770968). That repository represents a snapshot of the code and data at the time of publication. However, as needed, we will update the GitHub repository. 

In general, if you want to use the code for analysis, we recommend downloading the Zenodo repository, as it will work out of the box. If major changes are made to the code, we will update the GitHub repo and mirror them by incrementing the version in the Zenodo repository.

### Data compressed on GitHub
Some files are provided in a compressed format, but NOT in their uncompressed format, to facilitate easy access to the underlying data. For completeness, those files are listed here, with the format `<original filename>.zip` in all cases. These are listed below:

* In `figure_1/data`, the following compressed files are provided
	* `2025_pLDDT_scores_SHPRD.tsv.zip`
	* `2025_disorder_scores_SHPRD.tsv.zip`
	* `disorder_scores_SHPRD.tsv.zip`
	* `pLDDT_scores_SHPRD.tsv.zip`
	* `conservation_scores_SHPRD.tsv.zip`
	* `yeast_sequence_dataset.fasta.zip`
* In `extended_data_4/other_examples/data`, the following compressed files are provided
	*  `s_cerevisiae_clean.fasta.zip`


### Data excluded on GitHub
Some files are not shared on GitHub at all due to their large size, but are available in the Zenodo repository. These are listed below:

* In `figure_2/output_data/`, the following files are excluded from the GitHub Repo
	* `all_aligned_seqs.pkl`
* In `figure_5/evolution/figures/individual_intermaps/`, the following files are excluded from the GitHub Repo
	* 	`individual_intermaps.zip`
*  In `figure_5/proteome_wide_analysis/idr_msas`, the following directories are excluded from the GitHub Repo
	*  `metapredict_v1/`
	*  `metapredict_v3/`
*  In `figure_5/proteome_wide_analysis/generated_data`,  the following files are excluded from the GitHub Repo
	*  `ysn2idrs_v1_v1_disorder_0.8_50.pkl`
	*  `ysn2idrs_v3_v3_disorder_0.8_50.pkl`
*  In `figure_6/6B`,  the following files are excluded from the GitHub Repo:
	* `fig6B.zip`  
* In `extended_data_5/ED5_A`,  the following files are excluded from the GitHub Repo:
	* `extended_data_5_A.zip`
* In `extended_data_5/ED5_B`,  the following files are excluded from the GitHub Repo:
	* `extended_data_5_B.zip`
	* `read1_sgd_abf1wt_kpni_dm_sc_mhigh_15m.fastq.gz`
	* `read2_sgd_abf1wt_kpni_dm_sc_mhigh_15m.fastq.gz`
* In `extended_data_5/ED5_C`,  the following files are excluded from the GitHub Repo:
	* `extended_data_5_C.zip`

## Contact
For any questions or concerns, please contact co-corresponding authors [Alex Holehouse](https://www.holehouselab.com/) and [Philipp Korber](https://www.med.lmu.de/bmc/en/research/research-groups/korber-lab/).

## Main text figures

### Figure 1
* A - Schematic of conservation vs. disorder in a single protein (no data)
* B - `figure_1/disorder_vs_conservation_correlation.ipynb`
* C - `figure_1/disorder_vs_conservation_correlation.ipynb`
* D - `figure_1/disorder_vs_conservation_correlation.ipynb`
* E - Schematic of plasmid shuffling assay (no data)
* F - Data from dilution plates (see supplementary figures and supplementary tables)
* G - Schematic showing constructs (data from dilution plates)
* H - `figure_1/abf1_idr2.fasta` - amino acid sequence provided


### Figure 2
* A - Schematic of chemical conservation vs. sequence conservation
* B - `figures_2/composition_conservation.ipynb`
* C - `figures_2/alignment_schematic_figure_2c.ipynb`
* D - growth data from dilution plates (see supplementary figures and supplementary tables). Tree structure from Feng et al. *Scientific Reports* 2017. 
* E - growth data from dilution plates (see supplementary figures and supplementary tables) with protein function taken from UniProt and/or SGD entries on the proteins listed.
	
### Figure 3
* A - Schematic showing folded domain and IDR; growth data from dilution plates (see supplementary figures and supplementary tables).
* B - 	Schematic showing what sequence "shuffling" means; growth data from dilution plates (see supplementary figures and supplementary tables).
* C - Schematic showing where subregions are being shuffled (schematic to scale); growth data from dilution plates (see supplementary figures and supplementary tables).
* D - `figure_3/linear_analysis.ipynb`
* E - Schematic highlighting relative positions of the essential motif, Abf1-G4-like subregion, and $Gal4G4 subregion.
* F - Files in `figure_3/em_structural_prediction/`
* G - growth data from dilution plates (see supplementary figures and supplementary tables).
* H - Schematic (no data)
* I - growth data from dilution plates (see supplementary figures and supplementary tables).
* J - growth data from dilution plates (see supplementary figures and supplementary tables).
* K - growth data from dilution plates (see supplementary figures and supplementary tables).
* L - growth data from dilution plates (see supplementary figures and supplementary tables).
* M - growth data from dilution plates (see supplementary figures and supplementary tables).

### Figure 4
* A - Schematic of motif and context (top) and chemical vs. sequence specificity (bottom)
* B - `figure_4/coarse_grained_simulation_analysis.ipynb`
* C - growth data from dilution plates (see supplementary figures and supplementary tables).
* D - `figure_4/coarse_grained_simulation_analysis.ipynb`
* E - growth data from dilution plates (see supplementary figures and supplementary tables).
* F - growth data from dilution plates (see supplementary figures and supplementary tables).
* G - growth data from dilution plates (see supplementary figures and supplementary tables).
* H - growth data from dilution plates (see supplementary figures and supplementary tables).
* I - `figure_4/composition_analysis_fig_4.ipynb`
* J - growth data from dilution plates (see supplementary figures and supplementary tables).
* K - growth data from dilution plates (see supplementary figures and supplementary tables).

### Figure 5
* See readme in `figure_5/`
* A - Schematic of FINCHES-based analysis 
* B - Schematic of in silico evolution, with individual intermaps generated using `figure_5/evolution/in_silico_evo.ipynb`
* C - `figure_5/evolution/signal_vs_noise_final.ipynb`
* D - Schematic showing error-prone PCR workflow
* E - Schematic of expected behavior 
* F (top) - `figure_5/evolution/error_prone_pcr_analysis.ipynb`
* F (bottom) - `figure_5/evolution/error_prone_pcr_vs_in_silico_final.ipynb`

### Figure 6
* See readme in `figure_6/`
* A - Schematic of MNase-seq
* B - `figure_6/6B/MNase_coverage_Abf1_fig6B.xlsx`
* C - Schematic of Anchor Away (no data)
* D - Schematic of ODM-seq (no data)
* E - `figure_6/6E/6E_Kubik_responder_sites.tsv`
* F - `figure_6/6F/`
* G - `figure_6/6F/`
* H - `Schematic of overall model (no data)

## Extended data

### Extended Data 1
* A - `figure_1/disorder_vs_conservation_correlation.ipynb`
* B - `extended_data_1/fig_s1_analysis.ipynb`
* B inset - `extended_data_1/fig_s1_analysis.ipynb`
* C - `extended_data_1/fig_s1_analysis.ipynb`
* D - `extended_data_1/fig_s1_analysis.ipynb`
* E - `extended_data_1/fig_s1_analysis.ipynb`
* F - `extended_data_1/fig_s1_analysis.ipynb`
* G - 

### Extended Data 2
* A - `extended_data_2/chip_figure.ipynb`
* B - `figure_2/composition_conservation.ipynb`
* C - `figure_2/composition_conservation.ipynb`
* D - `figure_2/composition_conservation.ipynb`
* E - `figure_2/composition_conservation.ipynb`
* F - `figure_2/composition_conservation.ipynb`
* G - `figure_2/composition_conservation.ipynb`

### Extended Data 3
* A - `extended_data_3/ortholog_composition.ipynb`
* B - `figures_2 /compositional_analysis.ipynb`
* C - `extended_data_3/sequence_analysis_dead_seqs.ipynb`
* D - `extended_data_3/patterning_analysis.ipynb`
* E - `extended_data_3/patterning_analysis.ipynb`
* F - `extended_data_3/patterning_analysis.ipynb`
* G - `extended_data_3/sequence_alignment.docx` (alignment using [Needle](https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle)).

### Extended Data 4
* A - Schematic from simulation snapshot 
* B - Parameter matrix 
* C - `figure_4/coarse_grained_simulation_analysis.ipynb`
* D - `figure_4/coarse_grained_simulation_analysis.ipynb` 
* E - Amino acid sequences
* F - `extended_data_4/rad16_rad7_intermaps.ipynb`
* G - `extended_data_4/rad16_rad7_intermaps.ipynb`
* H - `extended_data_4/other_examples/signal_vs_noise_other_examples.ipynb`

### Extended Data 5
* See readme in `extended_data_5/`
* A - `extended_data_5/MNase_coverage_Abf1_ED5_A_B_C.xlsx`
* B - `extended_data_5/MNase_coverage_Abf1_ED5_A_B_C.xlsx`
* C - `extended_data_5/MNase_coverage_Abf1_ED5_A_B_C.xlsx`
* D - `extended_data_5/ED5_D/ED5_D_Kubik_responder_sites.tsv`
* E - `extended_data_5/ED5_E/ED5_E_Kubik_responder_sites.tsv`
* F - `extended_data_5/ED5_F/ED5_F_RNA_seq_boxplot_stat_v2.html`

### Extended Data 6
* See readme in `extended_data_6/`
* A - `extended_data_6/ED6_A/ED6_A_fixed_plus_one_grid_by_rep.tsv`
* B - `extended_data_6/ED6_B/ED6_B_spikein.tsv`
* C - `extended_data_6/ED6_C/ED6_C_k_resp_grid_by_rep.tsv`


### Extended Data 7
* See readme in `extended_data_7/`
* A - `extended_data_7/ED7_A_B_C/ED7_A_venn_diagram_invivo.tsv`
* B - `extended_data_7/ED7_A_B_C/ED7_B_venn_diagram_pwm.tsv`
* C - `extended_data_7/ED7_A_B_C/ED7_C_venn_diagram_pwm_vs_invivo.tsv`
* D - `extended_data_7/ED7_D_E_F_G/ED7_D_Hahn_bound_only.tsv`
* E - `extended_data_7/ED7_D_E_F_G/ED7_E_Hahn_resp_only.tsv`
* F - `extended_data_7/ED7_D_E_F_G/ED7_F_G_Hahn_bound_and_resp.tsv`
* G - `extended_data_7/ED7_D_E_F_G/ED7_F_G_Hahn_bound_and_resp.tsv`

### Extended Data 8
* See readme in `extended_data_8/`
* A - `extended_data_8/ED8_A/`
* B - `extended_data_8/ED8_B/`
* C - `extended_data_8/ED8_C/`
* D - `extended_data_8/ED8_D/`
* E - `extended_data_8/ED8_E/`


### Extended Data 9
* A - `extended_data_9/fus_sequences.fasta`
* B - Schematic of a phase diagram
* C - Schematic of solving free energy surface for Flory-Huggins to build temperature-dependent phase diagrams
* D - `extended_data_9/flory_huggins_phase_diagrams.ipynb`
* E - `extended_data_9/flory_huggins_phase_diagrams.ipynb`
* F - `extended_data_9/flory_huggins_phase_diagrams.ipynb`
* G - `extended_data_9/flory_huggins_phase_diagrams.ipynb`
* H - `extended_data_9/AF-Q13148-F1-model_v4.pdb` - PDB structure for TDP-43 used to highlight helical residues 321-342
* I - growth data from dilution plates (see supplementary figures and supplementary tables).

### Extended Data 10
* A - `extended_data_1/fig_s1_analysis.ipynb`
* B - `extended_data_1/fig_s1_analysis.ipynb`
* C - `extended_data_10/ED10_C/ED10_C_genome.tsv` 

## Supplementary information

### Supplementary figures 1 and 2
Growth data are provided in the `viablity_data/` directory, with two large PDFs that report original images for growth plates.

### Supplementary tables
For completeness, we also provide all supplementary tables in the `supplementary_tables/` directory (`supplementary_tables_2025_11_30.xlsx)`.

## License 
**DATA** is shared under CC BY 4.0 (Creative Commons Attribution 4.0 International) [link](https://creativecommons.org/licenses/by/4.0/deed.en). 

**CODE** is shared under the MIT License [link](https://mit-license.org/).

Any software tools used in this analysis are subject to their own licenses. 


## References
Feng, B., Lin, Y., Zhou, L., Guo, Y., Friedman, R., Xia, R., Hu, F., Liu, C. & Tang, J. Reconstructing yeasts phylogenies and ancestors from whole genome data. Sci. Rep. 7, 15209 (2017).
  
