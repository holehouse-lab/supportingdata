# Reproducing the computational figures from Perets *et al.* (2026)

#### Last updated 2026-08-28

This directory contains code, data and figures for the computational panels of:

> Perets, E. A., Spies, J. A., Cheong, J. H., Shi, L., Asamoto, D. K., Holehouse, A. S., Kim, J. E., Min, W., Neu, J. & Yan, E. C. Y. *Chain collapse, reduced dielectric, and water release drive protein phase separation.* **Biophysical Journal** 125 (2026). doi:[10.1016/j.bpj.2026.07.018](https://doi.org/10.1016/j.bpj.2026.07.018)

A copy of the paper, including the supporting information, is in [`perets_2025.pdf`](perets_2025.pdf).

## A proposed model for condensate formation

Biomolecular condensates are often (although not always...) together by weak, transient, heterogeneous interactions, which makes them challenging to study: most experimental approaches need labels, and labels perturb the very free energy landscape one is trying to measure. This work therefore uses a suite of *label-free* optical scattering and vibrational spectroscopies - UV resonance Raman, confocal Raman, dynamic light scattering, terahertz time-domain spectroscopy and stimulated Raman scattering imaging - alongside deep-learning-based ensemble prediction, to look inside condensates of the N-terminal domain of the RNA DEAD-box helicase 4 (DDX4<sup>[1-245]</sup>). 

Our work here proposes a model in which three things happen when DDX4<sup>[1-245]</sup> undergoes phase separation.:

1. The chain **collapses**, with the monomer $\langle R_h \rangle$ dropping from 4.2 nm to 1.6 nm on crossing the binodal. 
2. The **dielectric constant falls** inside the condensed phase, from 78 to 67 in H<sub>2</sub>O and from 69 to 49 in D<sub>2</sub>O. 
3.  **Water is released**, with the O-H stretching signal about 2.5% lower in the condensate core and about 4% lower at the interface. 

Because a lower dielectric strengthens the cation-pi and pi-pi interactions that drive condensation in the first place, these events plausibly form a positive-feedback loop that drives condensate growth. 

This is the proposed model, and while we do provide evidence supporting the observations for the model, it's worth being clear that the model here is an inferred model, as opposed to one where the order and details of these events are directly observed. 

## Scope of this repository

The Holehouse lab contributed the computational panels; the spectroscopy and light scattering were done in the Yan, Kim, Min and Neu labs. **As such, only the computational figures have code here.** The table below covers every figure in the paper so it is clear what is and is not reproducible from this directory.

| figure | content | code here |
|--------|---------|-----------|
| **[Fig S2](Fig_S2/)** | **coil-to-globule curve, LAMMPS + Mpipi-GG** | **yes** |
| **[Fig S3](Fig_S3/)** | **all-atom local conformational properties (CAMPARI)** | **yes** |
| **[Fig S4](Fig_S4/)** | **all-atom global dimensions, $R_g$/$R_h$ vs asphericity** | **yes** |
| **[Fig S5](Fig_S5/)** | **the same for STARLING-derived ensembles** | **yes** |
| **[Fig S6](Fig_S6/)** | **intramolecular scaling maps, wild-type and variants** | **yes** |
| **[Fig S7](Fig_S7/)** | **STARLING-derived $R_h$ across variants** | **yes** |
| **[Fig S8](Fig_S8/)** | **comparison of $R_h$ from four methods** | **yes** |
| **[Fig S9](Fig_S9/)** | **inferred $R_h$ for 144 proteins, with DDX4 overlaid** | **yes** |

## Layout

**Each figure directory is self-contained**, holding its own code, its own small data files and its own output, with a notebook that runs top to bottom using that directory as the working directory.

The one exception is bulk simulation trajectories, which are shared through [`data/`](data/README.md) at the repository root rather than duplicated per figure. That directory is **excluded from git** - the CAMPARI trajectories alone are ~750 MB - so a fresh clone will not have it. Notebooks that need it resolve it as `../data` by default, honour `PERETS_DATA_DIR` if you keep the files elsewhere, and print an explicit "MISSING" message naming what they expected rather than failing obscurely. `data/README.md` documents what belongs there and how it was generated; the trajectories are available on reasonable request. These data ARE provided on the full Zenodo repo though!

Each directory follows the same shape:

```
perets_2025/
├── data/              shared bulk trajectories - git-ignored, see data/README.md
├── Fig_S2/
│   ├── README.md      what the figure shows, how to run it, how it differs from the published panel
│   ├── figure_S2_*.ipynb
│   ├── data/          this figure's own small inputs, with a README documenting provenance
│   └── figures/       notebook output (PDF and PNG)
├── Fig_S3/
└── Fig_S4/
```

Every notebook README carries a "differences from the published figure" section. Where a number or a curve has shifted since publication - because of a package update, a corrected calculation, or restyling done in a vector graphics editor - it is documented there rather than quietly reconciled.

## Software

Analyses use [SOURSOP](https://soursop.readthedocs.io/) (2.0.5+) for simulation analysis, [STARLING](https://github.com/idptools/starling) for ensemble generation, [metapredict](https://github.com/idptools/metapredict) for disorder prediction, and standard `numpy`/`matplotlib`/`seaborn`. Individual figure READMEs list what each notebook actually needs; several need nothing beyond `numpy` and `matplotlib`.

## Data availability

Per the paper's data availability statement, all data will be shared upon reasonable request. Bulk trajectories and force field input files are not tracked in git; [`data/README.md`](data/README.md) documents the simulation provenance for the shared trajectories, and each figure's own `data/README.md` documents its local inputs and what is not included.
