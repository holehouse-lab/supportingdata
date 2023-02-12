## Additional Supporting Information
###### Last updated 2023-02-11

This repository contains additional information for the manuscript 

***SOURSOP: A Python package for the analysis of intrinsically disordered protein ensembles***

Jared Lalmansingh, Alex Keeley, Kiersten M. Ruff, Rohit V. Pappu, Alex S. Holehouse

For more information on the simulation data, please contact Alex or Jared directly!

## Important links:

#### [1. SOURSOP documentation](https://soursop.readthedocs.io/)
#### [2. SOURSOP GitHub repo](https://soursop.readthedocs.io/)

## Contents
### `/notebooks`
A collection of notebooks offering examples of analysis 

## Data used 
The SOURSOP paper re-analyzes a set of simulations performed by several different authors (see references below). 

Links to full IDP ensemble trajectories for the CAMPARI-generated  Monte Carlo simulation ensembles used in this paper are provided below. If these ensembles are used, PLEASE CITE the original papers (no need to cite this paper).

To access these trajectories [please click here](https://www.dropbox.com/sh/0cnrsodeatrsk9n/AADozaFr6T7VfEa1wfRgtPAQa?dl=0)

For the molecular dynamics simulations using [Amber99-disp](https://github.com/paulrobustelli/Force-Fields), ensembles were taken from Robustelli et al. 2019, and the full trajectories can be obtained from D.E. Shaw after agreeing to their license agreement (and as such, we cannot re-distribute them). However, Excluded Volume (EV) simulations are provided.

For CAMPARI simulations full simulations are provided

The list below reflects the mapping of different trajectories used to different papers.

* p27 - Das, Huang and Phillips PNAS 2016 [CAMPARI, ABSINTH, MC]
* ash1 - Martin & Holehouse et al. JACS 2016 [CAMPARI, ABSINTH, MC]
* notch - Sherry et al. PNAS 2017 [CAMPARI, ABSINTH, MC]
* actr - Robustelli et al. PNAS 2018 [DESMOND, amber99-disp, MD]
* asyn – Robustelli et al. PNAS 2018 [DESMOND, amber99-disp, MD]
* drkN – Robustelli et al. PNAS 2018 [DESMOND, amber99-disp, MD]
* ntail - Robustelli et al. PNAS 2018 [DESMOND, amber99-disp, MD]
* ntl9 – Peran & Holehouse et al. PNAS 2019 [CAMPARI, ABSINTH, MC]
* hnRNPA1 - Martin, Holehouse & Peran et al. Science 2020 [CAMPARI, ABSINTH, MC]
* p53\_1\_91 – Holehouse & Sukenik JCTC 2020 [CAMPARI, ABSINTH, MC]
* p53\_1\_91_S15E\_T18E\_S20E – This paper [CAMPARI, ABSINTH, MC]

## References

Das, R. K., Huang, Y., Phillips, A. H., Kriwacki, R. W., & Pappu, R. V. (2016). Cryptic sequence features within the disordered protein p27Kip1 regulate cell cycle signaling. Proceedings of the National Academy of Sciences of the United States of America, 113(20), 5616–5621.

Martin, E. W., Holehouse, A. S., Grace, C. R., Hughes, A., Pappu, R. V., & Mittag, T. (2016). Sequence Determinants of the Conformational Properties of an Intrinsically Disordered Protein Prior to and upon Multisite Phosphorylation. Journal of the American Chemical Society, 138(47), 15323–15335.

Sherry, K. P., Das, R. K., Pappu, R. V., & Barrick, D. (2017). Control of transcriptional activity by design of charge patterning in the intrinsically disordered RAM region of the Notch receptor. Proceedings of the National Academy of Sciences of the United States of America, 114(44), E9243–E9252.

Robustelli, P., Piana, S., & Shaw, D. E. (2018). Developing a molecular dynamics force field for both folded and disordered protein states. Proceedings of the National Academy of Sciences of the United States of America, 115(21), E4758–E4766.

Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, O., & Raleigh, D. P. (2019). Unfolded states under folding conditions accommodate sequence-specific conformational preferences with random coil-like dimensions. Proceedings of the National Academy of Sciences of the United States of America, 116(25), 12301–12310.

Martin, E. W., Holehouse, A. S., Peran, I., Farag, M., Incicco, J. J., Bremer, A., Grace, C. R., Soranno, A., Pappu, R. V., & Mittag, T. (2020). Valence and patterning of aromatic residues determine the phase behavior of prion-like domains. Science, 367(6478), 694–699.

Holehouse, A. S., & Sukenik, S. (2020). Controlling Structural Bias in Intrinsically Disordered Proteins Using Solution Space Scanning. Journal of Chemical Theory and Computation, 16(3), 1794–1805.