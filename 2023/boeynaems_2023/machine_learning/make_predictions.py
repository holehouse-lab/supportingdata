import os
import numpy as np
import pandas as pd
from protfasta import read_fasta
import joblib
from jax_unirep import get_reps
from jax_unirep.utils import load_params, get_weights_dir

# Constants
SIZE = 1900
PARAMS_PATH = get_weights_dir(paper_weights=SIZE)
MODEL_PATH = 'final_svm_model_1900.pkl'
FASTA_PATH = 'new_sequences.fasta'
PARAMS = load_params(PARAMS_PATH, SIZE)[1]

def load_fasta(filepath):
    """
    Load sequences from a FASTA file.

    Args:
        filepath (str): Path to the FASTA file.

    Returns:
        list: List of sequences.
    """
    seqs = read_fasta(filepath,invalid_sequence_action="ignore")
    headers, sequences = [], []
    for header, sequence in seqs.items():
        headers.append(header)
        sequences.append(sequence)
    return headers, sequences

def generate_embeddings(sequences, params, size):
    """
    Generate embeddings for a list of sequences using UniRep.

    Args:
        sequences (list): List of sequences.
        params: UniRep parameters.
        size (int): Size of the mLSTM model.

    Returns:
        np.ndarray: h_avg embeddings.
    """
    h_avg, _, _ = get_reps(sequences, params=params, mlstm_size=size)
    return h_avg

def load_model(filepath):
    """
    Load a trained model from a file.

    Args:
        filepath (str): Path to the saved model file.

    Returns:
        model: Loaded model.
    """
    return joblib.load(filepath)

def predict(model, embeddings):
    """
    Make predictions using the loaded model and embeddings.

    Args:
        model: Trained model.
        embeddings (np.ndarray): Feature matrix.

    Returns:
        np.ndarray: Predictions.
    """
    return model.predict(embeddings)

def main():
    """
    Main function to load sequences, generate embeddings, load the model,
    and make predictions.
    """
    headers,sequences = load_fasta(FASTA_PATH)
    embeddings = generate_embeddings(sequences, PARAMS, SIZE)
    model = load_model(MODEL_PATH)
    predictions = predict(model, embeddings)

    for sequence, prediction in zip(sequences, predictions):
        print(f'Sequence: {sequence[:50]}... Prediction: {prediction}')

if __name__ == "__main__":
    main()

