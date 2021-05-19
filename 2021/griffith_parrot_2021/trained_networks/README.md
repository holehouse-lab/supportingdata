# Making predictions with PARROT: general usage

Included in this directory:
* 6 Pre-trained PARROT networks (*.pt* extension) for phosphosite, activation domain and nucleation propensity prediction
* 2 Example files with sequences on which we can make predictions using these trained networks

First, if you haven't already, install the PARROT commandline tool using pip:
	
	pip install idptools-parrot

The command use for making predictions with pre-trained networks is `parrot-predict`. The basic usage for this command is:

	parrot-predict <sequence_file> <trained_network> <output_file> -d <datatype> -c <num_classes> -nl <number_hidden_layers> -hs <hidden_vector_size> [--probabilistic-classification]

## `sequence_file`
This is a tab/whitespace-delimited file that lists the sequences you want to make predictions for. Each row specifies a single sequence in the format:
	
    <seq_ID>	<AA_sequence>

See `example_seqs.tsv` for reference.

## `trained_network`
This is the file location of your pre-trained PARROT network. By PyTorch convention, I use *.pt* to specify these files, but this is not required.

## `output_file`
This is the file name and location where your predictions will be saved to. Outputted predictions will be listed in the same order as the provided `sequence_file`, along with the corresponding sequence IDs and amino acid sequences for each.

## `-d datatype`
These remaining arguments require some knowledge of how your network was trained. "Datatype" refers to whether your network is predicting *per-sequence* or *per-residue* values. Accordingly, the `-d` flag must be followed by either "sequence" or "residues", without the quotes.

## `-c num_classes`
This argument specifies whether the network was designed for a *classification* or *regression* task. If regression (i.e. the network outputs continuous, real-number values), "num_classes" should be 1. If classification (i.e. the network outputs a discrete class label), "num_classes" should be the number of possible classes a sequence/residue can be.

## `-nl number_hidden_layers` & `-hs hidden_vector_size`
These two arguments indicate the hyperparameters of the pre-trained PARROT network. These should be noted from when the network was trained, but if forgotten, it is possible to infer what these parameters should be based on the error messages that are produced when incorrect hyperparameters are supplied.

## `--probabilistic-classification`
This optional argument can only be used for sequence classification tasks with two possible classes. Its function is to change a binary class label output (0 or 1) into a real-number in [0,1] that represents confidence of the prediction, i.e., a value of 0.003 means that the network is fairly certain that the sequence is class 0, while 0.43 is rather uncertain.

There are a few other optional arguments that can be supplied. For information on these, see the full documentation at https://idptools-parrot.readthedocs.io/en/latest/

# Examples

## Predictions of phospho-tyrosines:

	parrot-predict phosphoY_seqs.tsv Y_phospho_network.pt phospho_predictions.tsv -d sequence -c 2 -nl 2 -hs 10

## Predictions of activation domains:

**Network trained on Erijman et al (ADpred) data:**

	parrot-predict AD_seqs.tsv AD_prediction_network.pt AD_predictions_Erijman.tsv -d sequence -c 2 -nl 2 -hs 10

**Network trained on Sanborn et al (PADDLE) data:**

	parrot-predict AD_seqs.tsv paddleTrain_AD_network.pt AD_predictions_Sanborn.tsv -d sequence -c 1 -nl 2 -hs 20

# Specific hyperparameters of provided networks
## Phosphorylation site predictions:
**S_phospho_network.pt**, **T_phospho_network.pt**, & **Y_phospho_network.pt**
* -d sequence
* -c 2
* -nl 2
* -hs 10
* While not required for the network to produce output, these networks are optimized to make predictions on sequences that are 19aa long, with the residue of interest (S, T, Y) in the central position.

## Activation domain function predictions:
**AD_prediction_network.pt**
* -d sequence
* -c 2
* -nl 2
* -hs 10
* Trained on data from Erijman et al. Best optimized to run on 30aa peptides.

**AD_prediction_network.pt**
* -d sequence
* -c 1
* -nl 2
* -hs 20
* Trained on data from Sanborn et al. Best optimized to run on 53aa peptides.

## AB42 nucleation score predictions:
**AB42_nucleation_predictor.pt**
* -d sequence
* -c 1
* -nl 2
* -hs 20
* Trained on data from Seuma et al. Best optimized to run on 42aa amyloid beta sequences, but may work on non-AB42 peptides.
