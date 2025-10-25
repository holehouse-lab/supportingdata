## SAXS back calculation
The experimental SAXS data here comes from the work of Holla *et al.* [1]. Specifically, both STARLING and experimental data are shared in the following directory:

	dCh_minus/
	dTRBP/
	sNh_plus/
	sNrich/
	
Synthetic scattering data is calculated from the ensembles using FOXS [2]. The resulting average scattering curves in each of these subdirectories are in `average_curve.dat`. There are additional files in here:

* `average_curve_error.dat` - includes an analytical estimation of the q-dependent error
* `average_curve_noise.dat` - includes stochastic noise from the analytical q-dependent error, and the q-dependent error itself
* `logfile_autoSCRT.txt` - logfile 

Worth noting, these files (calculated back from the STARLING ensembles) are ALSO found in the `../ensembles/<sequence_name>` directories. 

In ADDITION to these files, here we have the **experimental** data for the labelled protein. For each sequence, we have a couple of files:

* `<sequence_name>.dat` (e.g. `dTRBP.dat`) - the experimental, background-subtracted SAXS scattering data (in the [q, I(q), I(q) error], format).
* `original_filename.txt` - the original filename from which the experimental scattering data was copied from the Holla et al. data sets
* `mff_rg_<sequence_name>.csv` (e.g. `mff_rg_dTRBP.csv`) - contain $R_g$ values inferred from the experimental data using the Molecular Form Factor (MFF) developed by Riback et al. [3]. This number can be reproduced by uploading the experimental scattering curve to the [MFF webserver](https://sosnick.uchicago.edu/SAXSonIDPs/).
		
		
		
## References
[1] Holla, A., Martin, E. W., Dannenhoffer-Lafage, T., Ruff, K. M., König, S. L. B., Nüesch, M. F., Chowdhury, A., Louis, J. M., Soranno, A., Nettels, D., Pappu, R. V., Best, R. B., Mittag, T. & Schuler, B. **Identifying sequence effects on chain dimensions of disordered proteins by integrating experiments and simulations**. JACS Au 4, 4729–4743 (2024). [LINK](https://doi.org/10.1021/jacsau.4c00673)
  

[2] Schneidman-Duhovny, D., Hammel, M., Tainer, J. A. & Sali, A. **Accurate SAXS profile computation and its assessment by contrast variation experiments**. Biophys. J. 105, 962–974 (2013). [LINK](https://doi.org/10.1016/j.bpj.2013.07.020)

[3] Riback, J. A., Bowman, M. A., Zmyslowski, A. M., Knoverek, C. R., Jumper, J. M., Hinshaw, J. R., Kaye, E. B., Freed, K. F., Clark, P. L. & Sosnick, T. R. **Innovative scattering analysis shows that hydrophobic disordered proteins are expanded in water**. Science 358, 238–241 (2017). [LINK](doi.org/10.1073/pnas.1609579114)
  

  
		