## Setup
##### Last updated 2022-07-01

To run this analysis in a completely new conda environment the following set of commands should get you where you need to go.

> **NB:** If you're not familiar with conda and its role in Python environment configuration, we recommend [reading up on this](https://carpentries-incubator.github.io/introduction-to-conda-for-data-scientists/) first before continuing.

For the full functionality you will need to use Python 3.8. This is because the package [sparrow](https://github.com/holehouse-lab/sparrow) currently depends on [parrot](https://github.com/idptools/parrot) which itself depends on GPy and PyTorch which have some funky version dependencies, although we hope to fix these in the coming weeks...

	# create a new conda environment. We recommend 3.8 because as of right now
	# there are some versioning issues with PARROT/torch/gpy which impact
	# sparrow which is needed to compite 
	
	conda create -n tf_analysis python=3.8
	
	conda activate tf_analysis
	
	pip install cython
	conda install numpy scipy matplotlib jupyter 
	pip install shephard
	pip install metapredict
	pip install git+https://github.com/holehouse-lab/sparrow.git
	
Once these are installed, from within the same environment you will be able to open the notebooks `prepare_data.ipynb` and `compute_idr_params.ipynb`	
	
## Gene name changes

For a few of the genes the names provided aren't the canonical gene names, so I had to change some of the names to also include those genes names. Specifically

	CCDC101 ==> CCDC101/SGF29
	KDELC2 ==> KDELC2/POGLUT3
	SDCCAG3 ==> SDCCAG3/ENTR1