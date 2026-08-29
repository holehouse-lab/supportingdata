# STARLING ensembles

This directory holds the STARLING-predicted conformational ensembles behind Figures S5, S6, S7 and S8. **The ensemble files themselves are not tracked by git** (they total ~3.1 GB); this README is, so that the contents and the commands that produced them stay on the record. The files are available on reasonable request.

## What STARLING is

STARLING is a deep-learning generative model that predicts conformational ensembles of intrinsically disordered regions directly from sequence. It produces an ensemble in minutes on a laptop, where the equivalent all-atom simulation (Figures S3 and S4) takes about a month per system on a cluster. The model generates *inter-residue distance maps* with a denoising diffusion probabilistic model (DDPM) operating in the latent space of a variational autoencoder, then optionally reconstructs three-dimensional coordinates from those distance maps by multidimensional scaling.

Ensembles for this paper were generated with **STARLING version 1.0.0**, obtained April 2025, using the VAE/U-Net implementation with DDPM weights:

```
model-kernel-epoch.47-epoch_val_loss.0.03.ckpt
model-kernel-epoch.99-epoch_val_loss.1.72.ckpt
```

### A note on STARLING implementations

This work was originally done in the spring of 2025 using STARLING V1, [available here](https://github.com/idptools/starling/releases/tag/v1.0.0). As of summer 2026 (this paper's publication date) the default STARLING version has has changed to V2, which is a different model. To regenerate these figures the V1 of STARLING is required, ALTHOUGH the results here are robust to the STARLING version (with small difference).

## The command

STARLING installs a command-line tool, `starling`. The general form used here was:

```bash
starling <sequence> -b 400 -c <n_conformers>
```

The flags that matter:

| flag | meaning |
|------|---------|
| `-c`, `--conformations` | number of conformers to generate (default 400) |
| `-b`, `--batch_size` | how many conformers are generated in parallel (default 100) |
| `-r`, `--return_structures` | also reconstruct and write 3D coordinates, not just distance maps |
| `-o`, `--output_directory` | where to write output (default: current directory) |
| `--outname` | filename prefix for a single-sequence run |
| `-s`, `--steps` | number of DDPM denoising steps (default 30) |
| `--ionic_strength` | ionic strength in mM for the prediction (default 150) |
| `-d`, `--device` | compute device; auto-detected by default |

## The files

| file | sequence | conformers | used by |
|------|----------|-----------|---------|
| `ddx4_wt_6000.starling` | DDX4<sup>[1-245]</sup> wild-type | 6,000 | Figure S5, Figure S8 |
| `ddx4_wt_6000.xtc` / `.pdb` | the same ensemble, as coordinates | 6,000 | Figure S5 (asphericity) |
| `ddx4_wt.starling` | DDX4<sup>[1-245]</sup> wild-type | 1,600 | Figure S6, Figure S7 |
| `ddx4_24RtoK.starling` | 24RtoK variant, DDX4<sup>[1-236]</sup> | 1,600 | Figure S6, Figure S7 |
| `ddx4_14FtoA.starling` | 14FtoA variant | 1,600 | Figure S6, Figure S7 |
| `ddx4_37EDtoA.starling` | 37EDtoA variant | 1,600 | Figure S6, Figure S7 |
| `sp_repeat_245.starling` | (SP)<sub>122</sub>S null model, 245 residues | 800 | Figure S6 |
| `sp_repeat_236.starling` | (SP)<sub>118</sub> null model, 236 residues | 800 | Figure S6 |

The conformer counts above are read back from the files themselves. Only the wild-type 6,000-conformer run has its exact invocation recorded in the paper's methods; for the others the `-c` value is inferred from the file contents, and everything else about the invocation was the same.

To regenerate, taking the wild-type 6,000-conformer ensemble as the example:

```bash
starling MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSRFEDGDSSGFWRESSNDCEDNPTRNRGFSKRGGYRDGNNSEASGPYRRGGRGSFRGCRGGFGLGSPNNDLDPDECMQRTGGLFGSRRPVLSGTGNGDTSQSRSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGESSDTQGPKVT \
    -b 400 -c 6000 -r --outname ddx4_wt_6000
```

The `-r` flag is what produces the `.xtc` and `.pdb` alongside the `.starling` file. Figure S5 needs those coordinates, because asphericity requires the gyration tensor and so cannot be recovered from a distance map. Everything else these notebooks compute - radius of gyration, hydrodynamic radius in either estimator, mean distance maps - comes from the distance maps alone and needs no structures.

The same command with `-c 1600` and the appropriate sequence produces each of the 1,600-conformer files, and with `-c 800` the two null models.

## The sequences

### Wild-type DDX4<sup>[1-245]</sup> (245 residues)

```
MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSRFEDGDSSGFWRESSNDCEDNPTRNRGFSKRGGYRDGNNSEASGPYRRGGRGSFRGCRGGFGLGSPNNDLDPDECMQRTGGLFGSRRPVLSGTGNGDTSQSRSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGESSDTQGPKVT
```

Note this is the native 245-residue domain. The experimental construct carries an additional N-terminal `GAMGS` left over from TEV cleavage (see the supplemental sequence figure in the paper); the simulations use the domain sequence alone.

### 24RtoK (236 residues)

All 24 arginines in the DDX4<sup>[1-236]</sup> construct replaced by lysine, disrupting cation-pi interactions. This variant uses the shorter 236-residue construct to match previously reported experimental work, which is why it needs its own length-matched null model.

```
MGDEDWEAEINPHMSSYVPIFEKDKYSGENGDNFNKTPASSSEMDDGPSKKDHFMKSGFASGKNFGNKDAGECNKKDNTSTMGGFGVGKSFGNKGFSNSKFEDGDSSGFWKESSNDCEDNPTKNKGFSKKGGYKDGNNSEASGPYKKGGKGSFKGCKGGFGLGSPNNDLDPDECMQKTGGLFGSKKPVLSGTGNGDTSQSKSGSGSEKGGYKGLNEEVITGSGKNSWKSEAEGGES
```

### 14FtoA (245 residues)

All 14 phenylalanines replaced by alanine, disrupting both cation-pi and pi-pi interactions.

```
MGDEDWEAEINPHMSSYVPIAEKDRYSGENGDNANRTPASSSEMDDGPSRRDHAMKSGAASGRNAGNRDAGECNKRDNTSTMGGAGVGKSAGNRGASNSRAEDGDSSGAWRESSNDCEDNPTRNRGASKRGGYRDGNNSEASGPYRRGGRGSARGCRGGAGLGSPNNDLDPDECMQRTGGLAGSRRPVLSGTGNGDTSQSRSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGESSDTQGPKVT
```

### 37EDtoA (245 residues)

All 37 aspartate and glutamate residues replaced by alanine, disrupting electrostatic interactions.

```
MGAAAWAAAINPHMSSYVPIFAKARYSGANGANFNRTPASSSAMAAGPSRRAHFMKSGFASGRNFGNRAAGACNKRANTSTMGGFGVGKSFGNRGFSNSRFAAGASSGFWRASSNACAANPTRNRGFSKRGGYRAGNNSAASGPYRRGGRGSFRGCRGGFGLGSPNNALAPAACMQRTGGLFGSRRPVLSGTGNGATSQSRSGSGSARGGYKGLNAAVITGSGKNSWKSAAAGGASSATQGPKVT
```

### SP-repeat null models

Serine-proline dipeptide repeats, used to normalise the scaling maps in Figure S6. An SP repeat is an expanded, uniformly interacting IDR with no sequence-specific structure, so it is a null model for a chain of that length with nothing interesting going on. A null model has to be the same length as the sequence it normalises, which is why there are two.

The 245-residue version is `(SP)`x122 followed by a final `S`; the 236-residue version is `(SP)`x118. They can be generated in Python with `'SP' * 122 + 'S'` and `'SP' * 118`.

## File format and what can be computed from it

A `.starling` file stores the ensemble as a stack of inter-residue distance maps - shape `(n_conformers, n_residues, n_residues)`, float32 - together with the reconstructed 3D structures when `-r` was used. The distance maps dominate the file size; the wild-type 6,000-conformer ensemble is 6,000 x 245 x 245 float32 values, which is why that one file is 1.5 GB.

Load an ensemble with:

```python
from starling import load_ensemble
ensemble = load_ensemble('ddx4_wt_6000.starling')
```

Most quantities come from the distance maps and so return essentially instantly:

```python
ensemble.radius_of_gyration()                  # per conformer
ensemble.hydrodynamic_radius()                 # Nygaard et al. estimator (default)
ensemble.hydrodynamic_radius(mode='kr')        # Kirkwood-Riseman estimator
ensemble.distance_maps(return_mean=True)       # mean distance map
```

Passing `ignore_structures=True` to `load_ensemble` skips parsing the structures, which speeds up loading when only these quantities are needed.

**Asphericity is the exception.** It requires the gyration tensor and therefore real coordinates, so Figure S5 reads `ddx4_wt_6000.xtc`/`.pdb` with SOURSOP instead. The ensemble object does expose a `.trajectory` attribute returning a SOURSOP `SSProtein`, but calling the gyration-tensor code on that topology currently fails with an `mdtraj` attribute mismatch, so the trajectory files are read directly.

## References

Novak, B. J. M., Lotthammer, J. M., Emenecker, R. J. & Holehouse, A. S. *Accurate predictions of intrinsically disordered protein ensembles with STARLING.* **Nature** 652, 240-250 (2026). doi:[10.1038/s41586-026-10141-2](https://doi.org/10.1038/s41586-026-10141-2)

Nygaard, M., Kragelund, B. B., Papaleo, E. & Lindorff-Larsen, K. *An efficient method for estimating the hydrodynamic radius of disordered protein conformations.* **Biophys. J.** 113, 550-557 (2017). doi:[10.1016/j.bpj.2017.06.042](https://doi.org/10.1016/j.bpj.2017.06.042)

Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S. *SOURSOP: a Python package for the analysis of simulations of intrinsically disordered proteins.* **J. Chem. Theory Comput.** 19, 5609-5620 (2023). doi:[10.1021/acs.jctc.3c00190](https://doi.org/10.1021/acs.jctc.3c00190)
