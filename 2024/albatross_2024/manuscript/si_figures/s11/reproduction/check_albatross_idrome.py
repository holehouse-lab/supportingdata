"""Compare the Rg predictions of ALBATROSS against the published values from the IDRome dataset.

Attempt to reproduce Figure S15 from the ALBATROSS paper.
"""

import pandas as pd
from sklearn.metrics import r2_score
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from sparrow import Protein
from sparrow.predictors import batch_predict

TARGET = "scaled_rg"
MODEL_VERSION = 2
DATASET_VERSION = "initial"  # "initial" or "current"


def _plot_comparison(true_values: np.ndarray, sampled_values: np.ndarray, title: str) -> tuple[plt.Figure, plt.Axes]:
    """Plot the comparison of true and sampled values.

    Args:
        true_values : True values to compare.
        sampled_values : Sampled values to compare.
        title : Title of the plot.
    """
    # Compute the min and max values.
    all_values = np.concatenate([true_values, sampled_values])
    min_val, max_val = np.min(all_values), np.max(all_values)

    rmse = np.sqrt(np.mean((true_values - sampled_values) ** 2))
    r2 = r2_score(true_values, sampled_values)

    fig, ax = plt.subplots()
    sns.scatterplot(x=true_values, y=sampled_values, ax=ax)
    text_kwargs = dict(transform=ax.transAxes, verticalalignment="top", horizontalalignment="left")
    ax.text(0.05, 0.95, f"RMSE: {rmse:.5f}\nR2: {r2:.5f}", **text_kwargs)
    ax.plot([min_val, max_val], [min_val, max_val], "k--")
    ax.grid()
    ax.set_title(title)
    ax.set_xlabel("IDRome Rg (A)")
    ax.set_ylabel("ALBATROSS Rg (A)")

    return fig, ax


def check_albatross() -> None:
    """Compare Rg predictions of ALBATROSS against the published values from the IDRome dataset."""
    # Load the IDRome dataset.
    #df_idrome = pd.read_csv(f"data/IDRome_DB_{DATASET_VERSION}.csv")
    df_idrome = pd.read_csv(f"IDRome_DB.csv")
    sequence_names, sequences = df_idrome["seq_name"].tolist(), df_idrome["fasta"].tolist()
    save_path = f"./IDRome_{DATASET_VERSION}_{TARGET}_V{MODEL_VERSION}.png"

    # Extract true Rg values (and convert to Angstroms).
    true_rg_values = np.array(df_idrome["Rg/nm"].tolist()) * 10

    # Convert the sequences to Protein objects and predict the Rg values.
    protein_object_dict = {sequence_id: Protein(sequence) for sequence_id, sequence in zip(sequence_names, sequences)}
    output_dict = batch_predict.batch_predict(protein_object_dict, network=TARGET, version=MODEL_VERSION)
    pred_rg_values = np.array([output_dict[sequence_name][1] for sequence_name in sequence_names])

    # Plot the comparison.
    fig, _ = _plot_comparison(true_rg_values, pred_rg_values, title=f"Dataset: {DATASET_VERSION} | Model: {TARGET} (V{MODEL_VERSION})")
    print(f"Saving figure to {save_path}")
    fig.savefig(save_path)


if __name__ == "__main__":
    check_albatross()
