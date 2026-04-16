# Alston et al. (2026) — Analysis Code and Data

##### Last updated 2026-04-16

This repository contains the analysis code, data, and notebooks for reproducing the figures in the manuscript :

Alston, J. J., Soranno, A. & Holehouse, A. S. Conserved molecular recognition by an intrinsically disordered region in the absence of sequence conservation. *bioRxiv* (2026). doi:10.1101/2023.08.06.552128  



## Repository Structure

```
NTD_RNA_2026/
├── Figure_1/    Structural overview of coronavirus nucleocapsid NTD–RBD constructs
├── Figure_2/    Coarse-grained simulations showing an inert disordered region suppresses RNA binding
├── Figure_3/    Charge clustering in the SARS-CoV-2 NTD determines binding enhancement
├── Figure_4/    Orthologous NTD–RBD proteins show conserved bound-state ensembles
└── Figure_Supp/ Supplementary analyses (individual binding fits, SASA, scramble analyses)
```

### Figure 1
Structural and sequence conservation analysis of the NTD and RBD across six coronavirus orthologs (SARS-CoV-2, MERS, OC43, HKU1, 229E, MHV1). Includes AlphaFold-predicted structures and surface charge analysis via ChimeraX. Data only (no notebooks).

### Figure 2
Coarse-grained molecular dynamics analysis (Mpipi forcefield) of NTD–RBD + RNA binding.

| Notebook | Description |
|----------|-------------|
| `Figure_2BC.ipynb` | Calculates center-of-mass distances between NTD–RBD and RNA; fits bimodal distributions to determine bound/unbound states |
| `Figure_2E.ipynb` | Generates relative binding affinity (K_A*) box plots with statistical comparisons across GS10, GS25, SCO2, and RBD constructs |

### Figure 3
Analysis of how NTD charge patterning and sequence order control RNA-binding affinity.

| Notebook | Description |
|----------|-------------|
| `Plotting_Tile.ipynb` | Plots binding affinity for tiling mutants (T1–T26) that systematically reposition residues 30–50 |
| `Plotting_Charge.ipynb` | Correlates scrambled NTD charge positioning with binding affinity using IWD+ clustering metrics |
| `IWDShuffleExample.ipynb` | Generates 100,000 shuffled NTD sequence variants to assess charge patterning |

### Figure 4
Bound-state ensemble analysis across six coronavirus orthologs using inter-residue distance scaling maps.

| Notebook | Description |
|----------|-------------|
| `G.Summed_Interaction_*.ipynb` (×6) | Generates distance scaling maps (bound vs. unbound) for each ortholog (229E, HKU1, MERS, MHV1, OC43, SCO2) |
| `G.Summed_Interaction_SCO2_SCO2_for pdbannotation.ipynb` | Creates VMD-compatible coloring files for residue-level RNA contact visualization |
| `sequence_ensemble_analysis/4A_ntd_ensembles_.ipynb` | Analyzes radius of gyration distributions of NTD ensembles across orthologs |

### Supplementary Figures

| Notebook | Description |
|----------|-------------|
| `Fitted_Binding_Individual.ipynb` | Plots individual binding affinity measurements for NTD–RBD constructs across orthologs |
| `SCO2_RBD_SASA_analysis.ipynb` | Solvent-accessible surface area (SASA) analysis for SCO2 RBD structures |
| `additional_scramble_analysis/plot_scrambles_ordered_in_matrix.ipynb` | Visualizes 172 scrambled NTD sequences as color-coded matrices ordered by binding affinity |
| `scrambles_2025/process_sequences.ipynb` | Extracts sequences from simulation input files into FASTA format |
| `scrambles_2025/scramble_analysis.ipynb` | Mann-Whitney statistical tests comparing binding affinities of region-frozen NTD variants |

## Dependencies

### Python Packages

The notebooks require Python 3 and the following packages:

| Package | Purpose |
|---------|---------|
| `numpy` | Numerical computing |
| `pandas` | Data analysis and manipulation |
| `matplotlib` | Plotting and visualization |
| `seaborn` | Statistical data visualization |
| `scipy` | Scientific computing (curve fitting, statistical tests) |
| `mdtraj` | Molecular dynamics trajectory analysis |
| `soursop` | Simulation analysis (SSTrajectory) |
| `sparrow` | Protein sequence analysis (charge, structural properties) |
| `protfasta` | FASTA file parsing |
| `metapredict` | Protein disorder prediction |
| `statannotations` | Statistical annotations on plots |

Install with pip:

```bash
pip install numpy pandas matplotlib seaborn scipy mdtraj soursop sparrow protfasta metapredict statannotations
```

### External Software

Some figures were generated or processed with the following tools (not required to run the notebooks):

- **AlphaFold** — Structure prediction (Figure 1)
- **ChimeraX** — Electrostatic surface visualization (Figure 1)
- **VMD** — Molecular visualization and rendering (Figure 4)
- **LAMMPS** — Molecular dynamics simulations (trajectory generation)

## Citation

If you use code or data from this repository, please cite Alston et al. (2026).

