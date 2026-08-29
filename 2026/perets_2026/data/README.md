# Shared simulation data

This directory holds the large simulation trajectories that several figure notebooks read. **It is excluded from git** (see `../.gitignore`), because the trajectories run to hundreds of megabytes each. Everything else in `perets_2025` is tracked normally, and every figure directory keeps its own small data files locally.

If you have cloned this repository from GitHub, this directory will be empty or absent. The trajectories are available on reasonable request - see the paper's data availability statement - and the notebooks print a clear "MISSING" message telling you what they expected to find.

## Expected layout

```
data/
├── campari/
│   ├── ddx4_ntd/     full.xtc, full.pdb    DDX4[1-245], 25,000 conformers  (~331 MB)
│   └── ev/           full.xtc, full.pdb    self-avoiding random coil, 28,837 conformers  (~385 MB)
└── starling/
    ├── ddx4_wt_6000.starling               DDX4[1-245] wild-type, 6,000 conformers  (~1.5 GB)
    ├── ddx4_wt_6000.xtc, .pdb              the same ensemble as coordinates  (~7 MB)
    ├── ddx4_wt.starling                    DDX4[1-245] wild-type, 1,600 conformers  (~389 MB)
    ├── ddx4_24RtoK.starling                24RtoK variant, DDX4[1-236], 1,600 conformers  (~361 MB)
    ├── ddx4_14FtoA.starling                14FtoA variant, 1,600 conformers  (~389 MB)
    ├── ddx4_37EDtoA.starling               37EDtoA variant, 1,600 conformers  (~389 MB)
    ├── sp_repeat_245.starling              245-residue SP-repeat null model, 800 conformers  (~194 MB)
    └── sp_repeat_236.starling              236-residue SP-repeat null model, 800 conformers  (~181 MB)
```

Notebooks resolve this location as `../data` by default. Set the environment variable `PERETS_DATA_DIR` to point somewhere else if you keep the trajectories on another disk.

| directory | used by |
|-----------|---------|
| `campari/ddx4_ntd/` | Figure S3 (panels B, C), Figure S4 (both panels) |
| `campari/ev/` | Figure S3 (panels B, C) |
| `starling/ddx4_wt_6000.*` | Figure S5 (both panels), Figure S8 |
| `starling/ddx4_wt.starling` and the three variants | Figure S6, Figure S7 |
| `starling/sp_repeat_*.starling` | Figure S6 |

## Where these trajectories came from

Both are all-atom Monte Carlo simulations of the DDX4 N-terminal domain, run with the [CAMPARI](http://campari.sourceforge.net/) simulation engine (version 2) in conjunction with the ABSINTH implicit solvent model (parameter set `abs_3.2_opls.prm`) and the solution ion parameters of Mao *et al.* Simulations used the standard CAMPARI movesets and Hamiltonian parameters, as previously reported. As given in the paper's *Materials and methods*:

* Simulations were carried out at 340 K to facilitate robust sampling.
* 5 mM NaCl, in a spherical droplet with a fixed radius of 249 Angstroms, chosen to minimise finite-size effects.
* 20 independent simulations were run, each consisting of 30,000,000 equilibration steps followed by 90,000,000 production steps, with conformations saved every 80,000 steps during production.
* The 20 simulations comprise 10 started from a fully helical conformation and 10 started from a random coil. Starting from both extremes is a sampling control: agreement between the two sets means the result does not depend on the starting conformation.

These conditions are somewhat distinct from those examined experimentally in this paper, and so the analysis focuses primarily on local backbone behaviour rather than long-range interactions.

`ddx4_ntd/full.xtc` is the pooled production ensemble across all 20 simulations. The topology (`full.pdb`) carries 247 residues: the 245-residue DDX4<sup>[1-245]</sup> construct plus the two terminal caps CAMPARI adds.

**A note on the conformer count.** The pooled trajectory distributed here contains 25,000 conformers, which is the number quoted in the published caption for Figure S4. The *Materials and methods* section of the paper quotes a final ensemble of 23,773 conformers. The figures were made from the 25,000-conformer trajectory provided here.

### The excluded-volume / self-avoiding random coil ensemble

`campari/ev/` holds the self-avoiding random coil (SARC) reference, also referred to as the excluded volume limit. These simulations use the identical all-atom representation and the same engine, but with all attractive interactions set to zero, so the chain cannot overlap itself but has no reason to prefer any conformation. It is therefore the most conformationally heterogeneous ensemble this sequence can adopt.

This is what makes Figure S3 panels B and C interpretable. Both are comparisons against the SARC ensemble, so a deviation reflects a genuine sequence-encoded conformational bias rather than a trivial consequence of chain connectivity. The SARC simulations were performed as described previously in Holehouse *et al.* (2015).


## The STARLING ensembles

STARLING is a deep-learning generative model that predicts conformational ensembles of disordered regions directly from sequence, in minutes on commodity hardware rather than the weeks of compute an all-atom simulation needs. Ensembles here were generated with STARLING version 1.0.0 (obtained April 2025) using the VAE/U-Net implementation with DDPM weights.

**[`starling/README.md`](starling/README.md) documents these files in full** - the exact `starling` command-line invocations and what each flag does, the sequence of every variant and null model, the file format, and what can and cannot be computed from a `.starling` file without reconstructed coordinates. In brief:

* There are **two wild-type ensembles**, and they are not interchangeable. `ddx4_wt_6000.starling` (6,000 conformers) is used for Figures S5 and S8, where the distributions themselves are shown. `ddx4_wt.starling` (1,600 conformers) is used for Figures S6 and S7, where it is compared against the three variants, which all have 1,600 conformers - so the comparison is like for like. The two agree closely, as they should: mean hydrodynamic radius 41.42 versus 41.33 Angstroms.
* The **three variants** each remove one class of intramolecular interaction: 24 arginines to lysine (cation-pi), 14 phenylalanines to alanine (cation-pi and pi-pi), and 37 acidic residues to alanine (electrostatics). All three are phase-separation deficient.
* The **SP-repeat null models** are serine-proline dipeptide repeats used to normalise the Figure S6 scaling maps. There are two lengths because the 24RtoK variant uses the shorter DDX4[1-236] construct, and a null model must match the length of the sequence it normalises.
* **Asphericity needs coordinates**, which is why `ddx4_wt_6000.xtc`/`.pdb` sit alongside the `.starling` file and why Figure S5 reads them with SOURSOP. Everything else the notebooks compute comes from the stored distance maps.

## Analysis

All analysis is done with [SOURSOP](https://soursop.readthedocs.io/) (2.0.5 or newer) and the [STARLING](https://github.com/idptools/starling) package. The notebooks read these files directly and compute everything they plot; no precomputed intermediates are used.

## References

Vitalis, A. & Pappu, R. V. *ABSINTH: a new continuum solvation model for simulations of polypeptides in aqueous solutions.* **J. Comput. Chem.** 30, 673-699 (2009). doi:[10.1002/jcc.21005](https://doi.org/10.1002/jcc.21005)

Vitalis, A. & Pappu, R. V. *Methods for Monte Carlo simulations of biomacromolecules.* **Annu. Rep. Comput. Chem.** 5, 49-76 (2009). doi:[10.1016/S1574-1400(09)00503-9](https://doi.org/10.1016/S1574-1400(09)00503-9)

Mao, A. H. & Pappu, R. V. *Crystal lattice properties fully determine short-range interaction parameters for alkali and halide ions.* **J. Chem. Phys.** 137, 064104 (2012). doi:[10.1063/1.4742068](https://doi.org/10.1063/1.4742068)

Holehouse, A. S., Garai, K., Lyle, N., Vitalis, A. & Pappu, R. V. *Quantitative assessments of the distinct contributions of polypeptide backbone amides versus side chain groups to chain expansion via chemical denaturation.* **J. Am. Chem. Soc.** 137, 2984-2995 (2015). doi:[10.1021/ja512062h](https://doi.org/10.1021/ja512062h)

Novak, B. J. M., Lotthammer, J. M., Emenecker, R. J. & Holehouse, A. S. *Accurate predictions of intrinsically disordered protein ensembles with STARLING.* **Nature** 652, 240-250 (2026). doi:[10.1038/s41586-026-10141-2](https://doi.org/10.1038/s41586-026-10141-2)

Nygaard, M., Kragelund, B. B., Papaleo, E. & Lindorff-Larsen, K. *An efficient method for estimating the hydrodynamic radius of disordered protein conformations.* **Biophys. J.** 113, 550-557 (2017). doi:[10.1016/j.bpj.2017.06.042](https://doi.org/10.1016/j.bpj.2017.06.042)

Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S. *SOURSOP: a Python package for the analysis of simulations of intrinsically disordered proteins.* **J. Chem. Theory Comput.** 19, 5609-5620 (2023). doi:[10.1021/acs.jctc.3c00190](https://doi.org/10.1021/acs.jctc.3c00190)
