"""

Script to perform batch correction (with/without covariates) using ComBat.

"""
# import libraries
from typing import List, Optional

import hdbscan
import pandas as pd
from combat.pycombat import pycombat

from src.data_analysis.utils import (
    compute_pca,
    plot_pca,
    get_batch_info_from_sample_names,
    remove_single_sample_batches,
)


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
        raise ValueError(
            "The number of samples should match the number of batch indices."
        )

    # check if every batch factor has at least 2 samples
    if any([batch.count(i) < 2 for i in set(batch)]):
        raise ValueError("Every batch factor should have at least 2 samples.")

    # If not None, check if the number of covariates match the number of samples
    if covs is not None:
        if len(df.columns) != len(covs):
            raise ValueError(
                "The number of samples should match the number of covariates."
            )

    df_co = pycombat(data=df, batch=batch, mod=covs, mean_only=False)
    return df_co


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


# 1. Import data and metadata
path = "/Users/enrique/projects/local/heart-proteomics"
df = pd.read_table(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed.tsv",
    index_col=0,
    header=0,
)
# import metadata
metadata = pd.read_table(
    f"{path}/data/meta_data/heart_samples_meta_data_with_covs.tsv",
    index_col=False,
    header=0,
)


# Plot PCA of uncorrected imputed data
# transpose the dataframe to get samples as rows and features as columns
df_pca = compute_pca(df.T, n_components=5)

# add batch information to the dataframe
df_pca["batch"] = df_pca.index.str.split("-").str[0]

# plot PC1 vs PC2 with batch information using seaborn
# put the legend outside the plot
# save the plot as a png file
plot_pca(
    df_pca,
    title="PCA plot of uncorrected data",
    output_file=f"1.pca_uncorrected.png",
)

# keep samples only in "LV", "RV", "LA" and "RA" tissues from metadata
heart_parts = ["LV", "RV", "LA", "RA"]
metadata = metadata[metadata["organism_part_abbrev"].isin(heart_parts)]
samples_to_keep = metadata["sample_id"].tolist()

# keep samples in df that are also in samples_to_keep
df = df[[s for s in df.columns if s in samples_to_keep]]

# 2. Apply batch correction with covariate information
# Before apply batch correction, filter out batches with just one sample (otherwise the batch correction will fail).
batch_index = get_batch_info_from_sample_names(df.columns.tolist())
df = remove_single_sample_batches(df, batch_index)

# get covariate information from metadata.
columns = df.columns.tolist()
metadata = metadata[metadata["sample_id"].isin(columns)]
# reorder metadata to match the order of columns in df
metadata = metadata.set_index("sample_id").reindex(columns, axis=0).reset_index()
# get the covariates from metadata as list
covariates_index = metadata["covariate_index"].tolist()

# apply batch correction
df_corrected = apply_batch_correction(df, batch_index, covs=covariates_index)

# plot PCA of corrected data
# transpose the dataframe to get samples as rows and features as columns
df_pca = compute_pca(df_corrected.T, n_components=5)
# add batch information to the dataframe
df_pca["batch"] = df_pca.index.str.split("-").str[0]

# plot PC1 vs PC2 with batch information using seaborn
# put the legend outside the plot
# save the plot as a png file
plot_pca(df_pca, title="PCA plot of corrected data", output_file=f"2.pca_corrected.png")

# Export the dataframe with corrected batches
df_corrected.to_csv(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed-BatchCorrected.tsv",
    sep="\t",
)
