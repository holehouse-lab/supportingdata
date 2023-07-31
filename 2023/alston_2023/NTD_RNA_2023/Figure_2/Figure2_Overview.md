Figure 2. An Inert Disordered Region Can Suppress a Folded Domains RNA Binding Ability. 

A.
A snapshot of the bound state from a GS10-RBD + rU25 simulation trajectory. Visualization was done with Protein Imager(Tomasello et al. 2020). Simulations utilize the Mpipi forcefield (52). The model represents both amino acids and nucleotides as single beads with specific amino acid-amino acid and amino acid-nucleotide interactions. Folded domains are rigid and both disordered regions and nucleic acids are dynamic.

G. Tomasello, I. Armenia, G. Molla, The Protein Imager: a full-featured online molecular viewer interface with server-side HQ-rendering capabilities. Bioinformatics 36, 2909–2911 (2020).

B.
Example in Jupyter-Notebook Figure_2.ipnyb
The distances between the COM of the GS10-RBD and rU25 are plotted over the course of the simulation. A distance threshold (black line) is determined in C (see also Methods) and plotted to delineate the bound and unbound frames.


C.
Example in Jupyter-Notebook Figure_2.ipnyb
COM-COM distances from B are plotted as a histogram and show a bimodal distribution that correlate with the bound and unbound states of the protein. The distributions are fitted with dual Gaussians. A distance threshold, which separates bound and unbound frames, is determined by minimizing the overlap of the two populations.


D.
Schematic

E.
Example in Jupyter-Notebook Figure_2E.ipnyb uses all raw data stored as a CSV to generate significance plots
An apparent binding affinity (KA) is calculated by utilizing the fraction of bound and unbound frames and Eq. 1. This is then converted to a relative apparent binding affinity (KA*) by normalizing all values by dividing by the KA calculated from the SCO2 NTD-RBD + rU25 simulations. The boxes extend from the lower to upper quartile values of the data. The line within each box represent the median. The whiskers represent a 95% confidence interval. Significance is determined by a Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction. p-value annotation legend: (ns: 5.00e-02 < p <= 1.00e+00), (*: 1.00e-02 < p <= 5.00e-02), (**: 1.00e-03 < p <= 1.00e-02), (***: 1.00e-04 < p <= 1.00e-03), (****: p <= 1.00e-04)
