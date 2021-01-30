# CAMPARI simulation files
##### Last updated 2020-01-30

`run.key` contains the baseline keyfile that can be edited to run with different CAMPARI simulations (and includes information on which keywords should be altered).

Simulations require a `seq.in` [SEQFILE](http://campari.sourceforge.net/V3/inputfiles.html#FMCSC_SEQFILE), and in addition simulations with folded domains require a starting [PDB structure](http://campari.sourceforge.net/V3/keywords.html#PDBFILE).

In addition, a [PSWFILE](http://campari.sourceforge.net/V3/keywords.html#PSWFILE) is required to prevent sampling of backbone dihedrals in folded domains.

