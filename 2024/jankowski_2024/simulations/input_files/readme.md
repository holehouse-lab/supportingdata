## Simulation information for Jankowski et al.
###### Last updated 2024-03-02

This directory contains key information for running all-atom CAMPARI simulations of residues 754 to 803 (FFD region) in [*N. crassa* FRQ long isoform](https://www.uniprot.org/uniprot/P19970).

Simulations were run using the ABSINTH implicit solvent model [1].

The sequences simulated here are:

* WT: `PIAKEVMEPSGLGGVLPDDHFVMLVTTRRVVRPILQRQLSRSTTSEDTAE`
* AA: `PIAKEVMEPSGLGGVLPDDHFVMLVTTAAVVRPILQRQLSRSTTSEDTAE`
* HH: `PIAKEVMEPSGLGGVLPDDHFVMLVTTHHVVRPILQRQLSRSTTSEDTAE`
* HH_PRO:`PIAKEVMEPSGLGGVLPDDHFVMLVTTH*H*VVRPILQRQLSRSTTSEDTAE` - `H*` means protonated histidine.
* FCD\_WT: `TPDGLLPHHIVMTDKEKKKLVVRRLEQLFTGKISGRNMQRNQSMPSMDAP`
* FCD\_KKK_AAA: `TPDGLLPHHIVMTDKEAAALVVRRLEQLFTGKISGRNMQRNQSMPSMDAP`


Simulations were run using V2 of CAMPARI, although equivalent results should be achieved for other versions, and we'd recommend using V4, ALTHOUGH this almost certainly will necessitate some small changes to the keyfile 

* Normal simulations were run using the `full_keyfile.key`
* Excluded volume (EV) simulations were run using the `ev_keyfile.key`

In both cases, the `abs3.5_opls_unofficial.prm` parameters were used, which consist of the standard ABSINTH forcfield with the updated ion parameters of [Mao et al](https://aip.scitation.org/doi/10.1063/1.4742068) [2].

An example of an input sequence file (`wildtype_seq.in`) is provided. For mutants, changing the appropriate arginine residues to alanine (ALA), neutral histidine (HIE) or protonated histidine (HIP) can be done with a standard text editor. If the peptide net charge changes (i.e. ARG to ALA or ARG to HIE) then a compensatory change in the number of free ions should also be made (i.e. addition of two NA+ ions).


## References
[1] Vitalis, A. & Pappu, R. V. ABSINTH: A new continuum solvation model for simulations of polypeptides in aqueous solutions. J. Comput. Chem. 30, 673–699 (2009).

[2] Mao, A. H., & Pappu, R. V. (2012). Crystal lattice properties fully determine short-range interaction parameters for alkali and halide ions. The Journal of Chemical Physics, 137(6), 064104.