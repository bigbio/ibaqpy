# import libraries
import logging
import os
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# from inmoose.pycombat import pycombat_norm
from sklearn.cluster._hdbscan import hdbscan
from sklearn.decomposition import PCA

from ibaqpy.ibaq.ibaqpy_commons import IBAQ_NORMALIZED, SAMPLE_ID, PROTEIN_NAME

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def folder_retrieval(folder: str) -> dict:
    """
    Retrieve SDRF and ibaq results from a folder.
    :param folder: Folder to retrieve SDRF and ibaq results.
    :return:
    """

    folder = folder + os.sep if not folder.endswith(os.sep) else folder
    results = {"sdrf": [], "ibaq": []}
    items = os.listdir(folder)
    for item in items:
        try:
            results["sdrf"].extend(
                [
                    f"{folder}{item}/{i}"
                    for i in os.listdir(f"{folder}{item}/")
                    if i.endswith(".sdrf.tsv")
                ]
            )
            results["ibaq"].extend(
                [
                    f"{folder}{item}/{i}"
                    for i in os.listdir(f"{folder}{item}/")
                    if i.endswith("ibaq.csv") or i.endswith("ibaq.parquet")
                ]
            )
        except Exception as e:
            logger.warning(f"Error: {e}")
            if item.endswith(".sdrf.tsv"):
                results["sdrf"].append(folder + item)
            elif item.endswith("ibaq.csv"):
                results["ibaq"].append(folder + item)
            else:
                pass
    if len(results["sdrf"]) == 0:
        raise SystemExit("No SDRF founded!")
    if len(results["ibaq"]) == 0:
        raise SystemExit("No ibaq results founded!")
    if len(results["sdrf"]) != len(results["ibaq"]):
        raise SystemExit("Number of SDRFs should be equal to ibaq results!")

    return results


def generate_meta(sdrf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate ibaqpy metadata from SDRF. Each metadata contains four columns:

    - sample_id: Sample ID from every dataset (source name).
    - batch: PXD of every dataset (source name).
    - tissue: Tissue name of tissue-based dataset (characteristics[organism part]).
    - tissue_part: Tissue part of tissue-based dataset (characteristics[organism part]).

    param sdrf_df: _description_
    return: pd.DataFrame
    """
    sdrf_df.columns = [col.lower() for col in sdrf_df.columns]
    pxd = sdrf_df["source name"].values[0].split("-")[0]
    organism_part = [
        col for col in sdrf_df.columns if col.startswith("characteristics[organism part]")
    ]
    if len(organism_part) > 2:
        print(
            f"{pxd} Please provide a maximum of 2 characteristics[organism part], one for tissue name and the other for tissue part!"
        )
        exit(1)
    elif len(organism_part) == 0:
        print("Missing characteristics[organism part], please check your SDRF!")
        exit(1)

    meta_df = sdrf_df[["source name"] + organism_part]
    meta_df = meta_df.drop_duplicates()

    if len(meta_df.columns) == 2:
        meta_df["tissue_part"] = None
        meta_df.columns = ["sample_id", "tissue", "tissue_part"]
    else:
        if sdrf_df[organism_part[0]].nunique() > sdrf_df[organism_part[1]].nunique():
            a = "tissue_part"
            b = "tissue"
        else:
            a = "tissue"
            b = "tissue_part"
        meta_df.rename(
            columns={
                "source name": "sample_id",
                organism_part[0]: a,
                organism_part[1]: b,
            },
            inplace=True,
        )

    meta_df["batch"] = pxd
    meta_df = meta_df[["sample_id", "batch", "tissue", "tissue_part"]]
    meta_df = meta_df.drop_duplicates()

    return meta_df


def fill_samples(df, proteins):
    """
    Fill missing samples with 0 for all proteins.
    :param df: dataframe with samples in columns and proteins in rows.
    :param proteins: proteins to be used as index.
    :return:
    """
    df = pd.pivot_table(df, index=PROTEIN_NAME, columns=SAMPLE_ID, values=[IBAQ_NORMALIZED])
    df = df.reindex(proteins)
    df.columns = [pair[1] for pair in df.columns]
    df.index.rename(None, inplace=True)
    return df


def split_df_by_column(df: pd.DataFrame, cov_index_col: str) -> List[pd.DataFrame]:
    """
    Split a DataFrame by unique values of a specified column.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame to be split.
    cov_index_col : str
        The name of the column to split the DataFrame by.

    Returns
    -------
    List[pd.DataFrame]
        A list of pandas DataFrames, each containing rows with the same value in the specified column.
    """
    # Check if cov_index_col is in df
    if cov_index_col not in df.columns:
        raise ValueError(f"'{cov_index_col}' is not a column in the provided DataFrame.")

    # Use list comprehension to create the list of dataframes
    df_split = [df_group for _, df_group in df.groupby(cov_index_col)]

    return df_split


def filter_missing_value_by_group(df_input, col, non_missing_percent_to_keep):
    """
    Filters the dataframe by keeping columns with at least a specified percent of non-missing values
    in each group.

    Parameters:
    df_input (pandas.DataFrame): The input dataframe.
    col (str): The name of the column to group by.
    non_missing_percent_to_keep (float): The minimum percentage of non-missing values to keep a column.

    Returns:
    pandas.DataFrame: The filtered dataframe.
    """
    return df_input.groupby(col, as_index=False).filter(
        lambda x: len(x) < non_missing_percent_to_keep * len(df_input)
    )


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
    samples = [s.split("-")[0] for s in sample_list]
    batches = list(set(samples))
    index = {i: batches.index(i) for i in batches}

    return [index[i] for i in samples]


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


# define a function for batch correction
def apply_batch_correction(
    df: pd.DataFrame, batch: List[int], covs: Optional[List[int]] = None
) -> pd.DataFrame:
    """
    Get a dataframe and a list of batch indices as input and
    return a batch corrected dataframe with pycombat.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with the data to apply batch correction. Expected to have samples as columns and features as rows.
    batch : list
    A list of batch indices.
    covs : list
    A list of covariates to be used for batch correction.

    Warning
    -------
    Every batch should have at least 2 samples.

    Returns
    -------
    df_corrected : pd.DataFrame
    A batch-corrected dataframe.

    """

    # check if the number of samples match the number of batch indices
    if len(df.columns) != len(batch):
        raise ValueError("The number of samples should match the number of batch indices.")

    # check if every batch factor has at least 2 samples
    if any([batch.count(i) < 2 for i in set(batch)]):
        raise ValueError("Every batch factor should have at least 2 samples.")

    # If not None, check if the number of covariates match the number of samples
    if covs:
        if len(df.columns) != len(covs):
            raise ValueError("The number of samples should match the number of covariates.")

    from inmoose.pycombat import pycombat_norm
    df_co = pycombat_norm(counts=df, batch=batch, covar_mod=covs, mean_only=False)
    return df_co


# function to compute clusters
def find_clusters(df, min_cluster_size, min_samples) -> pd.DataFrame:
    """
    Compute clusters for a given dataframe.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with the data to be batched corrected.
    min_cluster_size : int
    The minimum size of clusters.
    min_samples : int
    The minimum number of samples in a neighborhood for a point to be considered as a core point.

    Returns
    -------
    df_clusters : pd.DataFrame
    A dataframe with the cluster assignments.
    """

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="euclidean",
        cluster_selection_method="eom",
        allow_single_cluster=True,
        cluster_selection_epsilon=0.01,
    )
    clusterer.fit(df)
    df["cluster"] = clusterer.labels_

    return df


# Function to run the iterative outlier removal pipeline.
# This function applies sequentially the following steps:
# 1. Compute principal components
# 2. Find clusters using HDBSCAN
# 3. Remove outliers
# 4. Repeat steps 1-3 until no outliers are found
def iterative_outlier_removal(
    df: pd.DataFrame,
    batch: List[int],
    n_components: int = 5,
    min_cluster_size: int = 10,
    min_samples: int = 10,
    n_iter: int = 10,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Get a dataframe and a list of batch indices as input and
    return a batch corrected dataframe with pycombat.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with the data to be batch corrected.
    batch : list
    A list of batch indices.
    n_components : int
    Number of principal components to be computed.
    min_cluster_size : int
    The minimum size of clusters.
    min_samples : int
    The minimum number of samples in a neighborhood for a point to be considered as a core point.
    n_iter : int
    Number of iterations to be performed.
    verbose : bool
    Whether to print and plot the number of outliers for each iteration.

    Returns
    -------
    df_filtered : pd.DataFrame
    A filtered dataframe.
    """

    # repeat steps 1-3 until no outliers are found
    # or the maximum number of iterations is reached
    # print the number of outliers for each iteration
    # save a plot of the principal components for each iteration

    # make dict with columns as key and batch as value
    batch_dict = dict(zip(df.columns, batch))

    for i in range(n_iter):
        print("Running iteration: ", i + 1)

        # compute principal components
        # transpose the dataframe to get samples as rows and features as columns
        df_pca = compute_pca(df.T, n_components=n_components)

        # compute clusters
        df_clusters = find_clusters(
            df_pca, min_cluster_size=min_cluster_size, min_samples=min_samples
        )
        print(df_clusters)
        # remove outliers from original dataframe
        outliers = df_clusters[df_clusters["cluster"] == -1].index.tolist()
        df_filtered_outliers = df.drop(outliers, axis=1)
        print(f"Number of outliers in iteration {i + 1}: {len(outliers)}")
        print(f"Outliers in iteration {i + 1}: {str(outliers)}")

        # update batch_dict based on the filtered dataframe
        batch_dict = {col: batch_dict[col] for col in df_filtered_outliers.columns}

        df = df_filtered_outliers

        # plot principal components PC1 vs PC2
        # save the plot as a png file
        # print the number of outliers for each iteration in the plot
        if verbose:
            plot_pca(
                df_clusters,
                output_file=f"iterative_outlier_removal_{i + 1}.png",
                x_col="PC1",
                y_col="PC2",
                hue_col="cluster",
                title=f"Iteration {i + 1}: Number of outliers: {len(outliers)}",
            )

        # break the loop if no outliers are found
        if len(outliers) == 0:
            break

    return df
