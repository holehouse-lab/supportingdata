# Data for Figure S9

## `IDR_dimensions.csv`

A compilation of intrinsically disordered proteins and regions for which a radius of gyration has been measured experimentally by small-angle X-ray scattering (SAXS). 145 entries, one per row.

| column | contents |
|--------|----------|
| `ID` | identifier from the source compilation; values run 1-150 with gaps, so IDs are not row numbers |
| `Name` | protein or construct name as recorded in the source |
| `Uniprot` | UniProt accession, where one applies |
| `Rg (Å)` | experimentally measured radius of gyration, Angstroms |
| `Error (Å)` | reported uncertainty on that measurement; **blank for 10 entries** |
| `Seq` | amino acid sequence of the region the measurement describes |

The file is UTF-8 with a byte-order mark and the column headers contain the Angstrom character, so it must be read with `encoding='utf-8-sig'`.

The sequence column is load-bearing rather than descriptive. The Nygaard conversion from $R_g$ to $R_h$ depends on chain length, so each protein's length is taken from its own sequence; the sequences span 15 to 698 residues.

## Where these numbers came from

The compilation was assembled by Alston *et al.* (2023) as the reference dataset for evaluating the analytical Flory random coil model against experimental measurements of disordered proteins. The values are extant SAXS-derived radii of gyration drawn from the literature, not measurements made for this paper. This notebook uses that table unchanged and only applies the $R_g$-to-$R_h$ conversion.

### Notes on the contents

* **145 entries, and a published caption that says 144.** Seven pairs of entries share an identical sequence while being distinct measurements - phosphorylated and unphosphorylated forms of the same protein (Ash1/pAsh1, PAGE4_WT and its phosphorylated form, CTCF-R domain WT and phosphorylated), and constructs recorded under two names (Tau ht40, histatin, pertactin-NTD/PNt, redAFP/sfAFP). Depending on how those are counted the compilation contains 138 unique sequences, 145 rows, or the 144 stated in the caption. The notebook plots all 145 rows.
* **Ten blank error values.** These entries report an $R_g$ with no stated uncertainty. They are kept, since the error is not used in this figure - it is read as `nan` and reported, not dropped.

## The two DDX4 points

The blue and yellow points overlaid on the plot are not from this file. They are the DLS measurements reported in Figure 3 of the paper: $\langle R_h \rangle$ = 4.2 nm for monomeric DDX4<sup>[1-245]</sup> under non-phase-separating conditions, and 1.6 nm under phase-separating conditions, both at 245 residues. They are entered directly in the notebook.

## References

Alston, J. J., Ginell, G. M., Soranno, A. & Holehouse, A. S. *The analytical Flory random coil is a simple-to-use reference model for unfolded and disordered proteins.* **J. Phys. Chem. B** 127, 4746-4760 (2023). doi:[10.1021/acs.jpcb.3c01619](https://doi.org/10.1021/acs.jpcb.3c01619)

Nygaard, M., Kragelund, B. B., Papaleo, E. & Lindorff-Larsen, K. *An efficient method for estimating the hydrodynamic radius of disordered protein conformations.* **Biophys. J.** 113, 550-557 (2017). doi:[10.1016/j.bpj.2017.06.042](https://doi.org/10.1016/j.bpj.2017.06.042)
