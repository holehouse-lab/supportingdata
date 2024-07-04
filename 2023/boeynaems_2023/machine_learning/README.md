# To replicate the analysis please run the following commands:
conda create -n venom python=3.11 -y
conda activate venom

pip install jax-unirep seaborn matplotlib protfasta

python train_model.py

# To make new predictions with a custom fasta file, you can use the make_predictions.py script.
# you'll need to edit the FASTA_PATH in the make_predictions.py to match the same path and name as the script you want to process.
python make_predictions.py
