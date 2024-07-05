This repository contains code for predicting protein condensation using embeddings generated with UniRep and a Support Vector Machine (SVM) classifier.

## Installation

To replicate the analysis, please run the following commands:

```bash
conda create -n venom python=3.11 -y
conda activate venom

pip install jax-unirep seaborn matplotlib protfasta
```
# Training the Model
To train the model, run:
```python train_model.py```

# Making Predictions
To make new predictions with a custom FASTA file, you can use the make_predictions.py script.

Ensure your custom FASTA file is in the correct location.
Run the following command, specifying your FASTA file path, the trained model path, and the desired output path:
bash

```python make_predictions.py --fasta_path path_to_new_fasta.fasta --model_path final_svm_model_1900.pkl --output_path predictions.csv```

By default, the script uses the following paths:

new_sequences.fasta for the FASTA file.
final_svm_model_1900.pkl for the model file.
predictions.csv for the output file.

# Repository Structure
- train_model.py: Script to train the SVM model using embeddings generated with UniRep.
- make_predictions.py: Script to make predictions on new sequences using the trained SVM model.

# Requirements
Python 3.11
jax-unirep
seaborn
matplotlib
protfasta
joblib
pandas
numpy
scikit-learn
tqdm
