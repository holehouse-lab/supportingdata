# About

This directory contains code for in silico evolution and error-prone PCR analysis associated with figure 5 of the manuscript.

We note that the mutagenesis protocol used generates new protein sequences that maintain similar interaction profiles to a given input sequence when interacting with a fixed partner sequence. The protocol uses either the Mpipi or CALVADOS coarse-grained models to evaluate the interaction matrices of the sequences. The mutagenesis process involves generating libraries of mutants through random mutations at the nucleotide level (converting protein back to DNA and performing mutagenesis in DNA space before converting to and selecting those that closely resemble the wild type interaction profile based on a defined threshold. The module allows for multiple rounds of mutagenesis to accumulate mutations while preserving interaction characteristics.

We note we are preparing a manuscript describing this protocol in detail (Flynn  et al. in prep); if you feel this type of analysis is useful for your work, **please** reach out to us for our updated version of this code! As of right now the underlying mutagenesis code is built into our internal (private) computational stack. We can share this if you need it, but what is likely more useful is to share our new version (CACTUS-WREN).

* `in_silico_evo.ipynb` - code to actually perform in silico evolution, with selection for chemical specificity .

* `libraries/` - output from `in_silico_evo.ipynb` goes here, which contains FASTA-files of evolved sequences. Note we actually only use the `chemlib_v2_350.fasta` and `chemlib_650_null_28.fasta` files for the manuscript.

* `signal_vs_noise_final.ipynb` - code to analyze in silico evolved sequence libraries

* `error_prone_pcr_sequences/` - sequences obtained from error-prone PCR

* `error_prone_pcr_analysis.ipynb` - code for analysis of outcomes from eror-prone PRC

* `error_prone_pcr_vs_in_silico_final.ipynb` - code that plots in silico vs. in vitro mutagenesis
