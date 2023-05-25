"""
Script to perform iterative outlier removal combining PCA and HDBSCAN

"""
# import libraries
from typing import List

import hdbscan
import pandas as pd

from src.data_analysis.utils import (
    compute_pca,
    plot_pca,
    get_batch_info_from_sample_names,
)


# function to compute clusters
def find_clusters(df, min_cluster_size=10, min_samples=10) -> pd.DataFrame:
    """
    Compute clusters for a given dataframe.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with the data to be batched corrected.
    min_cluster_size : int
    The minimum size of clusters.
    min_samples : int
    The minimum number of sampl es in a neighborhood for a point to be considered as a core point.

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

        # remove outliers from original dataframe
        outliers = df_clusters[df_clusters["cluster"] == -1].index.tolist()
        df_filtered_outliers = df.drop(outliers, axis=1)
        print(f"Number of outliers in iteration {i + 1}: {len(outliers)}")

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


# 1. Import data and metadata
path = "/Users/enrique/projects/local/heart-proteomics"
df = pd.read_table(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed-BatchCorrected.tsv",
    index_col=0,
    header=0,
)
# import metadata
metadata = pd.read_table(
    f"{path}/data/meta_data/heart_samples_meta_data_with_covs.tsv",
    index_col=False,
    header=0,
)


# Apply iterative outlier removal on corrected data
# get batch indices from the columns names
batch_indexes = get_batch_info_from_sample_names(df.columns.tolist())
# apply iterative outlier removal
df_filtered_outliers = iterative_outlier_removal(
    df,
    batch_indexes,
    n_components=5,
    min_cluster_size=20,
    min_samples=20,
    n_iter=3,
)

# plot PCA of corrected data with outliers removed
# transpose the dataframe to get samples as rows and features as columns
df_pca = compute_pca(df_filtered_outliers.T, n_components=5)

# add batch information to the dataframe
df_pca["batch"] = df_pca.index.str.split("-").str[0]

# plot PC1 vs PC2 with batch information using seaborn
# put the legend outside the plot
# save the plot as a png file
plot_pca(
    df_pca,
    title="PCA plot of corrected data with outliers removed",
    output_file=f"3.pca_corrected_outliers_removed.png",
)

# Export the dataframe with outliers removed
df_filtered_outliers.to_csv(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed-BatchCorrected-OutlierRemoved.tsv",
    sep="\t",
)
