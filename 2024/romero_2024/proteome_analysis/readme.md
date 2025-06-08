# About

This readme describes the steps needed to conduct 

## Step 1: Build domains from proteome
We build domains in two different ways, with DODO and with CHAINSAW. We do this to ensure the conclusions we come to are not biased by the specific domain decomposition approach.

The code for domain construction can be found in `/build_domains`/. The notebooks here enable the construction of both the DODO domains and the CHAINSAW domains where – at this stage – the domains are defined as stand-alone PDB files in `data/domains_<x>` where the `<x>` can be one of `_chainsaw/`, `_chainsaw_extended/` or `_dodo/`.




## Step 2: Analyze with FINCHES
Next we will analyze the solvent accessible residues on the surfaces of each domain using FINCHES.



The input of this analysis is an individual PDB file for each domain (with whatever domain definition procedure was defined) and the output is a SHEPHARD Domains file 


#### Output file:
Written to: `proteome_analysis/data/shprd_files`

Actual files:

* `shprd_chainsaw_domains_0.4.tsv`
* 

\

### Where are files generated?

* `data/shprd_files/shprd_{mode_name}_per_res_SASA.tsv` - in `build_sasa_annotations/domain_per_res_SASA.ipynb`.