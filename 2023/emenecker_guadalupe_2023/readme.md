## Additional Supporting Information
###### Last updated 2023-10-15

## Overview
This GitHub directory contains information for the upcoming manuscript **Sequence-ensemble-function relationships for disordered proteins in live cells** by Ryan J. Emenecker<sup>\*</sup>, Karina Guadalupe<sup>\*</sup>, Nora M. Shamoon, Shahar Sukenik<sup>✉</sup>, Alex S. Holehouse<sup>✉</sup> (* these authors contributed equally).

## Contact
For any questions or concerns please contact [Alex](https://www.holehouselab.com/) or [Shahar](https://www.sukeniklab.com/) via electronic mail. 

## How to use GOOSE
GOOSE can be used in two ways:

1. [GOOSE Colab Notebook](https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing) 
2. Local installation via PyPI (installation below)

## Installation guide

First, setup [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html). Having gotten conda working, the following install instructions will 

	# create a new environment
	conda create -n goose  python=3.9
	
	# activate that environment
	conda activate goose
	
	# install key initial requirements (could also use pip)
	conda install cython numpy
	
	# install goose
	pip install git+https://github.com/idptools/goose.git
		
We strongly recommend using Python 3.9 although versions above or equal to 3.7 should work. GOOSE has been tested on Python 3.9, 3.10, and 3.11. 

For examples of how to use GOOSE [see our online documentation](https://goose.readthedocs.io/).

Typical installation time is 2 minutes.

### System requirements
No special hardware is required and GOOSE will run on both macOS or Linux. In principle GOOSE should also run on Windows although this has not been tested.

## Other links

* [GOOSE GitHub repository](https://github.com/idptools/goose) - open source code for GOOSE, along with installation instructions 
* [GOOSE Documentation](https://goose.readthedocs.io/) - Documentation on how to use GOOSE for the design of disordered proteins
* [GOOSE Colab Notebook](https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing) - Colab notebook enabling simple design problems to be solved in the cloud

## Content

`bioinformatics/` - scripts and data for bioinformatic analysis

##