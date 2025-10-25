## About

Ensembles here were generated using the [IDPSam colab notebook](https://colab.research.google.com/github/giacomo-janson/idpsam/blob/main/notebooks/idpsam_experiments.ipynb). The output of the notebook generates three files:

* `pept.ca.top.pdb` (topology file)
* `pept.ca.traj.dcd` (trajectory file)
* `pept.seq.fasta` (sequence file)

Ensembles here have 1000 conformers.

We then calculated the radius of gyration using the topology/trajectory file generating:

* `RG.csv` - instantaneous $R_g$ for every conformer.
* `RG_mean.csv` - average $R_g$ value.
* `RG_std.csv` - standard deviation of $R_g$ values.

We note that ensembles here were generated at coarse-grained resolution, many because all-atom resolution is recovered using the [cg2aa](https://academic.oup.com/bioinformatics/article/32/8/1235/1744650) program, which means any all-atom detail we'd see is really a product of cg2aa and not from IDPSam.
