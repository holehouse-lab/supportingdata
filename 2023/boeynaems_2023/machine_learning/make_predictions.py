import os
import numpy as np
import pandas as pd
from protfasta import read_fasta
import joblib
from jax_unirep import get_reps
from jax_unirep.utils import load_params, get_weights_dir
import argparse

# Constants
SIZE = 1900
PARAMS_PATH = get_weights_dir(paper_weights=SIZE)
PARAMS = load_params(PARAMS_PATH, SIZE)[1]

def load_fasta(filepath):
    """
    Load sequences from a FASTA file.

    Args:
        filepath (str): Path to the FASTA file.

    Returns:
        list: List of sequences.
    """
    seqs = read_fasta(filepath, invalid_sequence_action="ignore")
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

def save_predictions(headers, sequences, predictions, output_path):
    """
    Save the predictions to a CSV file.

    Args:
        headers (list): List of sequence headers.
        sequences (list): List of sequences.
        predictions (np.ndarray): Predictions from the model.
        output_path (str): Path to save the CSV file.
    """
    results = pd.DataFrame({
        'Header': headers,
        'Sequence': sequences,
        'Prediction': predictions
    })
    results.to_csv(output_path, index=False)

def main():
    """
    Main function to load sequences, generate embeddings, load the model,
    and make predictions.
    """
    parser = argparse.ArgumentParser(description="Make predictions on new sequences using a trained SVM model.")
    parser.add_argument('--fasta_path', type=str, default='new_sequences.fasta', help="Path to the FASTA file containing new sequences.")
    parser.add_argument('--model_path', type=str, default='final_svm_model_1900.pkl', help="Path to the trained SVM model file.")
    parser.add_argument('--output_path', type=str, default='predictions.csv', help="Path to save the predictions CSV file.")

    args = parser.parse_args()

    headers, sequences = load_fasta(args.fasta_path)
    embeddings = generate_embeddings(sequences, PARAMS, SIZE)
    model = load_model(args.model_path)
    predictions = predict(model, embeddings)

    save_predictions(headers, sequences, predictions, args.output_path)

if __name__ == "__main__":
    main()

