import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.svm import SVC
from sklearn.metrics import classification_report,confusion_matrix
import joblib
from jax_unirep import get_reps
from jax_unirep.utils import load_params, get_weights_dir

# Constants
SIZE = 1900
CLASSES = 2
DATA_PATH = "data/mass_spec_{}_class_labeled_dataset.tsv".format(CLASSES)
TOP_10_FILEPATH = "data/data.tsv"
PARAMS_PATH = get_weights_dir(paper_weights=SIZE)

def load_data(filepath, classes):
    """
    Load dataset from a specified file path.

    Args:
        filepath (str): Path to the dataset file.
        classes (int): Number of classes in the dataset.

    Returns:
        pd.DataFrame: Loaded dataset.
    """
    return pd.read_csv(filepath, sep="\t")

def load_model_params(filepath, size):
    """
    Load model parameters.

    Args:
        filepath (str): Path to the parameters directory.
        size (int): Size of the mLSTM model.

    Returns:
        params: Loaded parameters.
    """
    return load_params(filepath, size)[1]

def generate_embeddings(sequences, params, size):
    """
    Generate embeddings for a list of sequences using UniRep.

    Args:
        sequences (list): List of sequences.
        params: UniRep parameters.
        size (int): Size of the mLSTM model.

    Returns:
        np.ndarray, np.ndarray, np.ndarray: h_avg, h_final, c_final embeddings.
    """
    return get_reps(sequences, params=params, mlstm_size=size)

def load_or_generate_embeddings(data, params, size):
    """
    Load pre-generated embeddings if available, otherwise generate them.

    Args:
        data (pd.DataFrame): Data containing sequences.
        params: UniRep parameters.
        size (int): Size of the mLSTM model.

    Returns:
        np.ndarray, np.ndarray, np.ndarray: h_avg, h_final, c_final embeddings.
    """
    if os.path.exists("h_avg_{}.npy".format(size)):
        h_avg = np.load("h_avg_{}.npy".format(size))
        h_final = np.load("h_final_{}.npy".format(size))
        c_final = np.load("c_final_{}.npy".format(size))
    else:
        sequences = data['Sequence'].tolist()
        h_avg, h_final, c_final = generate_embeddings(sequences, params, size)
        np.save("h_avg_{}.npy".format(size), h_avg)
        np.save("h_final_{}.npy".format(size), h_final)
        np.save("c_final_{}.npy".format(size), c_final)
    return h_avg, h_final, c_final

def train_model(features, targets, param_grid, kfold):
    """
    Train an SVM model with cross-validation and grid search.

    Args:
        features (np.ndarray): Feature matrix.
        targets (pd.Series): Target vector.
        param_grid (dict): Parameter grid for grid search.
        kfold (StratifiedKFold): Stratified K-Fold cross-validator.

    Returns:
        GridSearchCV: Fitted GridSearchCV object.
    """
    svc_model = SVC(kernel='linear', gamma='auto', cache_size=500)
    grid_search = GridSearchCV(svc_model, param_grid, cv=kfold, return_train_score=True)
    grid_search.fit(features, targets)
    return grid_search

def evaluate_model(model, features, targets, kfold):
    """
    Evaluate the model using cross-validation.

    Args:
        model (SVC): Trained SVC model.
        features (np.ndarray): Feature matrix.
        targets (pd.Series): Target vector.
        kfold (StratifiedKFold): Stratified K-Fold cross-validator.

    Returns:
        np.ndarray, np.ndarray: True labels and predictions from all folds.
    """
    predictions = []
    true_labels = []
    for train_index, test_index in kfold.split(features, targets):
        X_train, X_test = features[train_index], features[test_index]
        y_train, y_test = targets[train_index], targets[test_index]
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        predictions.append(y_pred)
        true_labels.append(y_test)

    all_predictions = np.concatenate(predictions)
    all_true_labels = np.concatenate(true_labels)
    return all_true_labels, all_predictions

def generate_report(true_labels, predictions, targets):
    """
    Generate a classification report.

    Args:
        true_labels (np.ndarray): True labels from cross-validation.
        predictions (np.ndarray): Predictions from cross-validation.
        targets (pd.Series): Original target vector.
    """
    print("Classification Report:")
    print(classification_report(true_labels, predictions, labels=np.unique(targets)))

def predict_top_10(model, top_10_data, params, size):
    """
    Predict classes for top 10 sequences.

    Args:
        model (SVC): Trained SVC model.
        top_10_data (pd.DataFrame): Data containing top 10 sequences and their classes.
        params: UniRep parameters.
        size (int): Size of the mLSTM model.
    """
    top_10_havg, _, _ = get_reps(top_10_data["sequence"].tolist(), params=params, mlstm_size=size)
    top_10_predictions = model.predict(top_10_havg)
    for new, orig in zip(top_10_predictions, top_10_data["class"]):
        print(new, orig)

def plot_confusion_matrix(true_labels, predictions, targets):
    """
    Plot and save the confusion matrix.

    Args:
        true_labels (np.ndarray): True labels from cross-validation.
        predictions (np.ndarray): Predictions from cross-validation.
        targets (pd.Series): Original target vector.
    """
    cm = confusion_matrix(true_labels, predictions, labels=np.unique(targets))
    plt.figure(figsize=(10, 8))
    sns.set(style="whitegrid")
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False, 
                xticklabels=np.unique(targets), yticklabels=np.unique(targets),
                linewidths=.5)
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig('confusion_matrix.png', dpi=300)
    plt.show()

def save_model(model, filepath):
    """
    Save the trained model to a file.

    Args:
        model (SVC): Trained SVC model.
        filepath (str): Path to save the model.
    """
    joblib.dump(model, filepath)

def main():
    """
    Main function to load data, train and evaluate the model, and save the results.
    """
    data = load_data(DATA_PATH, CLASSES)
    params = load_model_params(PARAMS_PATH, SIZE)
    h_avg, _, _ = load_or_generate_embeddings(data, params, SIZE)
    targets = data['aggCluster={}'.format(CLASSES)]

    kfold = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)
    param_grid = {"C": [0.1 * i for i in range(1, 20)]}

    grid_search = train_model(h_avg, targets, param_grid, kfold)
    best_C = grid_search.best_params_['C']

    final_model = SVC(C=best_C, kernel='linear', gamma='auto', cache_size=500, probability=True)
    final_model.fit(h_avg, targets)

    all_true_labels, all_predictions = evaluate_model(final_model, h_avg, targets, kfold)
    generate_report(all_true_labels, all_predictions, targets)

    plot_confusion_matrix(all_true_labels, all_predictions, targets)


    top_10 = pd.read_csv(TOP_10_FILEPATH, sep="\t", names=["sequence", "class"])
    predict_top_10(final_model, top_10, params, SIZE)

    save_model(final_model, 'final_svm_model_1900.pkl')

if __name__ == "__main__":
    main()

