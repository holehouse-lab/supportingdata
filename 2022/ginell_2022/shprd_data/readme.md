# About these data (Human)

The data in this directory come from several sources.

## Human proteome
Two copies of the human proteome are provided here, and both are based on the reference human [proteome obtained from UniProt](https://www.uniprot.org/proteomes/UP000005640).


`human_proteome_validated.fasta` - contains all protein sequences with those that have non-canonical amino acids removed. This prevents exceptions with other tools later in various analysis pipelines. This contains 20,393 proteins.

In addition, there's a a SHEPHARD proteins file (`shprd_proteins_filtered_human_af2_f1acc.tsv`) which contains the protein sequence information for the proteins where there are AlphaFold2 annotations. This contains 20,061 proteins (i.e. ~300 fewer than the full dataset).

## Structural information
All structural annotations (discussed below) are based on predicted 3D structures from AlphaFold2.

### Structural context
Structural context for every residue was calculated using the DISICL protocol for residue structural classification (Nagy et al 2014). This defines the local phi/psi context into one of 19 distinct structural classes. DISICL values were calculated using code in the SESCA package (Nagy et al 2019).

These assignments can be found in `shprd_tracks_DISICL_human.tsv`

### Binding accessibility
Binding accessibility scores are based on SASA calculated for each residue using a probe radius of 7 A (in contrast to the 1.5 usually used for water). These values were calculated using [MDTraj](https://mdtraj.org/).

These scores can be found in `shprd_tracks_filtered_af2_f1_solv_acc.tsv`.

### DSSP assignment
DSSP assignments for each AlphaFold2 residue are provided in `shprd_tracks_dssp_human.tsv`. The DSSP assignments were calculated using using [MDTraj](https://mdtraj.org/), based on the original Kabsch and Sander algorithm.

These scores can be found in `shprd_tracks_dssp_human.tsv`.

#### References
Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., … Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583–589.

Tunyasuvunakool, K., Adler, J., Wu, Z., Green, T., Zielinski, M., Žídek, A., Bridgland, A., Cowie, A., Meyer, C., Laydon, A., Velankar, S., Kleywegt, G. J., Bateman, A., Evans, R., Pritzel, A., Figurnov, M., Ronneberger, O., Bates, R., Kohl, S. A. A., … Hassabis, D. (2021). Highly accurate protein structure prediction for the human proteome. Nature, 596(7873), 590–596.

Nagy, G., & Oostenbrink, C. (2014). Dihedral-based segment identification and classification of biopolymers I: proteins. Journal of Chemical Information and Modeling, 54(1), 266–277.

Nagy, G., Igaev, M., Jones, N. C., Hoffmann, S. V., & Grubmüller, H. (2019). SESCA: Predicting Circular Dichroism Spectra from Protein Molecular Structures. Journal of Chemical Theory and Computation, 15(9), 5087–5102.

McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M., Hernández, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J., & Pande, V. S. (2015). MDTraj: a modern, open library for the analysis of molecular dynamics trajectories. Biophysical Journal, 109(8), 1528–1532.

Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637.

## Disorder predictions
Disorder predictions found in `shprd_domains_idrs_metapredict_v2.tsv` come from [metapredict](https://metapredict.readthedocs.io/en/latest/) (V2). This is the dataset we recommend using, although we also include IDRs as predicted with metapredict V1 in `shprd_domains_human_IDRs_0.420.tsv` (0.42 here is the cutoff used for deliniating IDRs and non-IDRs, which is a conservative cutoff minimizing false positives but potentially lead to false negatives).

To make these predictions install metapredict (`pip install metapredict`) and then simply run:

	metapredict-predict-idrs human_proteome_validated.fasta  --mode shephard-domains-uniprot  -o shprd_domains_idrs_metapredict_v2.tsv --verbose
	
This will generate a SHEPHARD-compliant domains file as output, using the UniProt IDs as the unique IDs.	

#### References
Emenecker, R. J., Griffith, D., & Holehouse, A. S. (2021). Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure. Biophysical Journal, 120(20), 4312–4319.

Emenecker, R. J., Griffith, D., & Holehouse, A. S. (2022). Metapredict V2: An update to metapredict, a fast, accurate, and easy-to-use predictor of consensus disorder and structure. In bioRxiv (p. 2022.06.06.494887). https://doi.org/10.1101/2022.06.06.494887


	metapredict-predict-idrs human_proteome_validated.fasta  --mode shephard-domains-uniprot  -o shprd_domains_idrs_metapredict_v2.tsv --verbose


## Human polar-rich low-complexity sequences
Human QSGNTP-rich low complexity domains (also called polar-rich low-complexity domains) were generated using [sparrow](https://github.com/idptools/sparrow). LCDs were identified using the LCD-definition method developed by Guitierrez et al


These domains are found in `shprd_domains_human_QSGNTP_LCD.tsv`.

#### References
Gutierrez, J. I., Brittingham, G. P., Karadeniz, Y., Tran, K. D., Dutta, A., Holehouse, A. S., Peterson, C. L., & Holt, L. J. (2022). SWI/SNF senses carbon starvation with a pH-sensitive low-complexity sequence. eLife, 11. https://doi.org/10.7554/eLife.70344


## Prion-like domains

Prion-like domains (PLDs) were defined using default settings in the PLAAC algorithm.

These domains can be found here: `shprd_domains_human_plaac_PLDs.tsv`

### References

Lancaster, A. K., Nutter-Upham, A., Lindquist, S., & King, O. D. (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. Bioinformatics , 30(17), 2501–2502.

Alberti, S., Halfmann, R., King, O., Kapila, A., & Lindquist, S. (2009). A systematic survey identifies prions and illuminates sequence features of prionogenic proteins. Cell, 137(1), 146–158.

## Protein abundance data
Found in `shprd_protein_attributes_concentrations.tsv` comes from the supplementary information in:

#### References
Hein, M. Y., Hubner, N. C., Poser, I., Cox, J., Nagaraj, N., Toyoda, Y., Gak, I. A., Weisswange, I., Mansfeld, J., Buchholz, F., Hyman, A. A., & Mann, M. (2015). A human interactome in three quantitative dimensions organized by stoichiometries and abundances. Cell, 163(3), 712–723.


