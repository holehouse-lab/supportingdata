# IDRs from other predictors
###### Last updated 2021-05-31

**metapredict** contains a stand-alone function for extracting IDRs from disorder profiles.

While - by default - this function uses the metapredict disorder scores, we also include the functionality to extract IDRs from _other_ disorder profiles as well.

The tutorial below walks through how one could use this approach with ODiNPred to extract out IDRs. However, if you can convert your per-residue disorder scores into a Python list, any predictor could be used here

## IDRs from ODiNPred
[ODiNPred](https://st-protein.chem.au.dk/odinpred) is a fantastic disorder predictor that combines NMR data with evolutionary analysis to make data-driven predictions of disorder.

To use ODiNPred, shimmy on over to their web server add your sequences and your email address and hit **submit**. Once the prediction has run, ODiNPred will email you a per-protein file with the per-residue disorder scores named something like `DisorderPredictions....txt` with the following format.

	 Residue  No. Zscore Disorder-Probability
	 M         1  1.6758  0.9890
	 E         2  1.6679  0.9890
	 E         3  1.7146  0.9886
	 P         4  1.5467  0.9957
	 Q         5  1.4447  0.9961
	 S         6  1.6193  0.9954
	...		  ...   ...  	...
	
Say our file was named 	`DisorderPredictionssp_P04637_P53_HUMANC.txt` then the code in the notebook here (`idrs_from_odinpred.ipynb`) walks step-by-step through the process of extracting out the IDRs based on the ODiNPred scores.

Note that we have no affiliation with the ODiNPred folks - it's just a great disorder predictor! 