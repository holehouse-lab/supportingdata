Figure 1. Coronavirus nucleocapsid proteins possess a disordered, poorly-conserved N-terminal domain (NTD) and a more well-conserved folded RNA binding domain (RBD). 


A.
Schematic showing full-length nucleocapsid protein architectures from coronaviruses. The nucleocapsid protein contains three IDRs (NTD, Linker, CTD) and two folded 
domain (RBD, and Dimerization domains).





B.
Per-residue conservation calculated over 45 orthologous NTR-RBD constructs, including SCO2, MERS, OC43, HKU1, 229E and MHV1. Conservation is calculated based on the positional Shannon entropy, with values shown only for residues where 80% or more of orthologous posses a residue. The NTD contains many gaps in a relatively poor alignment, while the RBD is almost uniformly populated with relatively highly-conserved residues.




C.
Overlay of RBD structures for SCO2, MERS, OC43, HKU1, 229E and MHV1, revealing a high degree of structural conservation in the RBD fold. 
To generate structures 5 Alphafold structures were generated with Alphafold Colab (https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)
the Most similar structural conformation for each ortholog was chosen by assessing RMSD values between each conformation. The 5 structures chosen
are deposited in the 'Structure' folder.

SCO2-
MERS-
OC43-
HKU1-
229E-
MHV1-




D.
Surface charge properties of the six RBD structures overlaid in panel C, highlighting differences in surface charge properties despite conservation of overall fold. 
Charge was modeled using ChimeraX using Molecule_Display, Electrostatic, Coulombic with Amber 20 recommended default charges and atom types for standard residues

Coulombic values for OC43_5 (1).pdb_A SES surface #7.1: minimum, -15.74, mean -1.01, maximum 11.67
Coulombic values for MHV1_5 (1).pdb_A SES surface #8.1: minimum, -10.64, mean 0.50, maximum 12.18
Coulombic values for MERS_1 (1).pdb_A SES surface #9.1: minimum, -8.77, mean 1.55, maximum 10.84
Coulombic values for HKU1_5 (1).pdb_A SES surface #10.1: minimum, -10.56, mean 0.84, maximum 11.66
Coulombic values for 229E_3 (1).pdb_A SES surface #11.1: minimum, -13.31, mean 2.95, maximum 14.66
Coulombic values for SCO2_2 (1).pdb_A SES surface #12.1: minimum, -9.48, mean 2.01, maximum 15.15