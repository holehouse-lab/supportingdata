# Organisms in `aligned.fasta`

`aligned.fasta` contains 12 FRQ (Frequency clock protein) orthologs from filamentous fungi in the class Sordariomycetes, aligned to a common length of 1073 columns. The headers are locus tags from Ensembl Fungi annotations rather than species names; the `-t26_1`-style suffixes are transcript/isoform identifiers, not part of the gene name.

| Header | Organism | Ungapped length |
|---|---|---|
| `sp\|P19970\|FRQ_NEUCR` | *Neurospora crassa* OR-74A, UniProt reference entry | 989 |
| `NEUTE1DRAFT_67585-t26_1` | *Neurospora tetrasperma* | 989 |
| `NEUDI_101601T0` | *Neurospora discreta* | 989 |
| `SMAC_03705-t26_1` | *Sordaria macrospora* | 927 |
| `FGRAMPH1_01T19375` | *Fusarium graminearum* PH-1 | 993 |
| `FVEG_04686-t26_1` | *Fusarium verticillioides* | 992 |
| `FFUJ_13245-t31_1` | *Fusarium fujikuroi* | 960 |
| `FPRN_05785-t42_1` | *Fusarium proliferatum* ET1 | 992 |
| `FOC1_g10014445-t38_1` | *Fusarium oxysporum* f. sp. *cubense* race 1 | 991 |
| `TRIREDRAFT_121670-t26_1` | *Trichoderma reesei* | 1014 |
| `TRIVIDRAFT_20727-t45_1` | *Trichoderma virens* | 1007 |

The sampling is four Sordariaceae sequences (three *Neurospora*, one *Sordaria*), five *Fusarium*, and two *Trichoderma*, spanning the Hypocreales and Sordariales.

## Caveat for conservation calculations

*N. crassa* is present twice: `NCU02265-t26_1` and the UniProt `P19970` entry are the same protein (both 989 residues ungapped), so it is effectively double-weighted in any column-wise frequency calculation.
