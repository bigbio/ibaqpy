# import libraries
from typing import List, Optional

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from combat.pycombat import pycombat
from sklearn.decomposition import PCA


# function to compute principal components
def compute_pca(df, n_components=5) -> pd.DataFrame:
    """
    Compute principal components for a given dataframe.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with samples as rows and features as columns.
    n_components : int
    Number of principal components to be computed.

    Returns
    -------
    df_pca : pd.DataFrame
    A dataframe with the principal components.
    """

    pca = PCA(n_components=n_components)
    pca.fit(df)
    df_pca = pca.transform(df)

    df_pca = pd.DataFrame(
        df_pca, index=df.index, columns=[f"PC{i}" for i in range(1, n_components + 1)]
    )

    return df_pca


# get batch info from sample names
def get_batch_info_from_sample_names(sample_list: List[str]) -> List[int]:
    """
    Expected as input a list of sample names with SDRF-like format: {PRIDE_PROJECT_ID}-{SAMPLE_ID}
    and return a list of batch indices (a.k.a. factor levels)

    :param sample_list: list of sample names
    :return: list of batch indices
    """
    batches = [s.split("-")[0] for s in sample_list]
    factor_levels = np.array(batches)
    unique_levels, batch_indexes = np.unique(factor_levels, return_inverse=True)
    return batch_indexes.tolist()


# define a function to remove batches with only one sample.
# takes as input a dataframe with samples in columns and protein IDs in rows, and a list of batch indices
# returns a dataframe with batches with only one sample removed
def remove_single_sample_batches(df: pd.DataFrame, batch: list) -> pd.DataFrame:
    """
    Remove batches with only one sample.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with samples in columns and protein IDs in rows.
    batch : list
    A list of batch indices.

    Returns
    -------
    df_filtered : pd.DataFrame
    A filtered dataframe.
    """

    # make dict with columns as key and batch as value
    batch_dict = dict(zip(df.columns, batch))

    # from the batch_dict, get the batches with only one sample
    single_sample_batch = [
        k for k, v in batch_dict.items() if list(batch_dict.values()).count(v) == 1
    ]

    # remove batches with only one sample
    df_single_batches_removed = df.drop(single_sample_batch, axis=1)

    # update batch_dict based on the filtered dataframe
    # batch_dict = {col: batch_dict[col] for col in df_single_batches_removed.columns}

    return df_single_batches_removed


def plot_pca(
    df_pca,
    output_file,
    x_col="PC1",
    y_col="PC2",
    hue_col="batch",
    palette="Set2",
    title="PCA plot",
    figsize=(8, 6),
):
    """
    Plots a PCA scatter plot and saves it to a file.

    Args:
    df_pca (pd.DataFrame): DataFrame containing PCA data.
    output_file (str): File name to save the plot as an image.
    x_col (str, optional): Column name for x-axis. Defaults to "PC1".
    y_col (str, optional): Column name for y-axis. Defaults to "PC2".
    hue_col (str, optional): Column name for hue (grouping variable). Defaults to "batch".
    palette (str, optional): Color palette for the plot. Defaults to "Set2".
    title (str, optional): Title for the plot. Defaults to "PCA plot".
    figsize (tuple, optional): Figure size as (width, height) in inches. Defaults to (5, 5).
    """

    # Create a new figure with the specified size
    fig, ax = plt.subplots(figsize=figsize)

    # Create a scatterplot using seaborn
    sns.scatterplot(x=x_col, y=y_col, hue=hue_col, data=df_pca, palette=palette, ax=ax)

    # Set the plot title, x-axis label, and y-axis label
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)

    # Set the legend location and adjust the bounding box
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

    # Adjust the layout to fit the legend within the figure
    plt.tight_layout()

    # Save the plot as an image file
    plt.savefig(output_file, bbox_inches="tight")
