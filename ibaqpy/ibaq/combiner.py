import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ibaqpy.ibaq.ibaqpy_commons import load_feature, load_sdrf
from ibaqpy.ibaq.imputation_methods import impute_missing_values
from ibaqpy.ibaq.utils import (
    compute_pca,
    get_batch_info_from_sample_names,
    generate_meta,
    folder_retrieval,
    filter_missing_value_by_group,
    split_df_by_column,
    fill_samples,
    iterative_outlier_removal,
    plot_pca,
    remove_single_sample_batches,
    apply_batch_correction,
)

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


class Combiner:
    def __init__(self, data_folder: os.PathLike, covariate: str = None, organism: str = "HUMAN"):
        """
        Initialize a Combiner instance to process and combine SDRF and iBAQ data.

        Parameters:
        data_folder : os.PathLike
            Path to the folder containing SDRF and iBAQ files.
        covariate : str, optional
            Covariate to be used in data processing, default is None.
        organism : str, optional
            Organism filter for protein names, default is "HUMAN".

        Raises:
        FileNotFoundError
            If the specified data folder does not exist or is not a directory.

        Notes
        -----
        This method initializes various attributes for data processing, retrieves
        SDRF and iBAQ files from the specified folder, and filters protein data
        by the specified organism.
        """
        self.df_pca = compute_pca(self.df_corrected.T, n_components=5)
        self.df_corrected = None
        self.batch_index = get_batch_info_from_sample_names(self.df.columns.tolist())
        self.df_pca = None
        self.df_filtered_outliers = None
        self.batch_index = get_batch_info_from_sample_names(self.df.columns)
        self.samples_number = None
        self.datasets = None
        self.samples = self.df.columns.tolist()
        self.proteins = self.df["ProteinName"].unique().tolist()
        logger.info("Combining SDRFs and ibaq results ...")
        self.data_folder = Path(data_folder)
        if not self.data_folder.exists() or not self.data_folder.is_dir():
            raise FileNotFoundError(f"Data folder {self.data_folder} does not exsit!")
        self.covariate = covariate
        files = folder_retrieval(str(self.data_folder))
        self.metadata, self.df = pd.DataFrame(), pd.DataFrame()
        for sdrf in files["sdrf"]:
            sdrf_df = load_sdrf(sdrf)
            self.metadata = pd.concat([self.metadata, generate_meta(sdrf_df)])
        self.metadata = self.metadata.drop_duplicates()
        self.metadata.index = self.metadata["sample_id"]

        for ibaq in files["ibaq"]:
            self.df = pd.concat([self.df, load_feature(ibaq)])
        self.df = self.df[self.df["ProteinName"].str.endswith(organism)]
        self.df.index = self.df["SampleID"]
        self.df = self.df.join(self.metadata, how="left")
        logger.info(self.metadata, self.df.head)

    def read_data(self, meta: str, ibaq: str, organism="HUMAN", covariate=None):
        """
        Reads and processes iBAQ and metadata files, filtering protein data by organism.

        Parameters:
        meta : str
            Path to the metadata CSV file.
        ibaq : str
            Path to the iBAQ CSV file.
        organism : str, optional
            Organism filter for protein names, default is "HUMAN".
        covariate : str, optional
            Covariate to be used in data processing, default is None.

        Notes
        -----
        The method updates the instance's dataframe and metadata attributes by
        reading the specified files, filtering proteins by the given organism,
        and joining metadata to the iBAQ data.
        """

        self.covariate = covariate
        self.df = pd.read_csv(ibaq, index_col=0)
        self.metadata = pd.read_csv(meta)
        self.df = self.df[self.df["ProteinName"].str.endswith(organism)]
        self.df.index = self.df["SampleID"]
        self.metadata = self.metadata.drop_duplicates()
        self.df = self.df.join(self.metadata, how="left")

    def imputer(self, covariate_to_keep: list = None):
        """
        Impute missing values in the combined iBAQ results DataFrame.

        This method processes the DataFrame by filtering, filling, and imputing
        missing values based on specified covariates. It ensures that only columns
        with a sufficient percentage of non-missing values are retained and performs
        imputation using KNN or other specified methods.

        Parameters:
        covariate_to_keep : list, optional
            A list of covariate values to retain in the DataFrame. Only rows with
            these covariate values will be kept.

        Raises:
        SystemExit
            If the specified covariate contains fewer than two unique values.

        Notes
        -----
        - The method modifies the instance's DataFrame by imputing missing values
          and potentially altering its structure.
        - The imputation process requires samples as columns and proteins as rows.
        """
        logger.info("Imputing merged ibaq results ...")
        # Keep only columns 'sample_id' and covariate from df_metadata
        if self.covariate:
            if len(self.metadata[self.covariate].unique()) < 2:
                raise SystemExit(
                    f"{self.covariate} should contain at least two different covariates!"
                )

        # Keep only rows within covariate_to_keep, you can keep tissue or tissue part you want.
        if covariate_to_keep:
            self.df = self.df[self.df[self.covariate].isin(covariate_to_keep)]

        # keep columns with at least 30% of non-missing values in each covariate_index group
        self.df = filter_missing_value_by_group(
            self.df, col="ProteinName", non_missing_percent_to_keep=0.3
        )

        # TODO: Data for imputation should take samples as columns, proteins as rows. [Expression Matrix]
        # Also need to fill the proteins didn't show in original results for each sample.
        if self.covariate:
            # split df by covariates
            df_list = split_df_by_column(self.df, cov_index_col=self.covariate)
            df_list = [fill_samples(df, self.proteins) for df in df_list]

            # impute missing values with KNNImputer for every df in df_list
            df_list = impute_missing_values(df_list)

            # concatenate all dataframes in df_list into one dataframe
            self.df = pd.concat(df_list, axis=1)
        else:
            self.df = fill_samples(self.df, self.proteins)
            self.df = impute_missing_values(self.df)

        self.datasets = list(set([sample.split("-")[0] for sample in self.samples]))
        logger.info(f"DataFrame head after imputation:\n{self.df.head()}")

    def outlier_removal(
        self,
        n_components: int = None,
        min_cluster_size: int = None,
        min_samples_num: int = None,
        n_iter: int = None,
    ):
        """
        Remove outliers from the imputed data using an iterative approach and plot the PCA results.

        This method applies iterative outlier removal on the imputed data, updates the filtered
        DataFrame, and generates a PCA plot of the corrected data with outliers removed.

        Parameters:
        n_components : int, optional
            Number of principal components to compute. Defaults to a third of the unique batch indices.
        min_cluster_size : int, optional
            Minimum size of clusters for outlier detection. Defaults to the median number of samples per batch.
        min_samples_num : int, optional
            Minimum number of samples in a neighborhood for a point to be considered a core point.
            Defaults to the median number of samples per batch.
        n_iter : int, optional
            Number of iterations for outlier removal. Defaults to 5.

        Notes
        -----
        - The method modifies the instance's DataFrame by removing outliers.
        - A PCA plot is saved as 'pca_corrected_outliers_removed.png'.
        """
        logger.info("Removing outliers from imputed data ...")
        # Apply iterative outlier removal on imputed data
        # get batch indices from the columns names
        batches = [sample.split("-")[0] for sample in self.samples]
        self.samples_number = {dataset: batches.count(dataset) for dataset in self.datasets}
        min_samples = round(np.median(list(self.samples_number.values())))
        if min_samples == 1:
            min_samples = 2
        # apply iterative outlier removal
        self.df_filtered_outliers = iterative_outlier_removal(
            self.df,
            self.batch_index,
            n_components=(n_components if n_components else round(len(set(self.batch_index)) / 3)),
            min_cluster_size=min_cluster_size if min_cluster_size else min_samples,
            min_samples=min_samples_num if min_samples_num else min_samples,
            n_iter=n_iter if n_iter else 5,
        )
        logger.info(self.df_filtered_outliers)
        # plot PCA of corrected data with outliers removed
        # transpose the dataframe to get samples as rows and features as columns
        self.df_pca = compute_pca(
            self.df_filtered_outliers.T,
            n_components=(n_components if n_components else round(len(set(self.batch_index)) / 3)),
        )

        # add batch information to the dataframe
        self.df_pca["batch"] = self.df_pca.index.str.split("-").str[0]

        # plot PC1 vs PC2 with batch information using seaborn
        # put the legend outside the plot
        # save the plot as a png file
        plot_pca(
            self.df_pca,
            title="PCA plot of corrected data with outliers removed",
            output_file="pca_corrected_outliers_removed.png",
        )

    def batch_correction(self, n_components: int = None, tissue_parts_to_keep: int = None):
        """
        Apply batch effect correction to the data and plot PCA results.

        This method performs batch correction on the data using specified covariates
        and plots PCA before and after correction. It filters out batches with only
        one sample and optionally retains specific tissue parts.

        Parameters:
        n_components : int, optional
            Number of principal components to compute. Defaults to a third of the unique batch indices.
        tissue_parts_to_keep : int, optional
            Number of tissue parts to retain in the data.

        Notes
        -----
        - The method modifies the instance's DataFrame by applying batch correction.
        - PCA plots are saved as 'pca_uncorrected.png' and 'pca_corrected.png'.
        """
        logger.info("Applying batch effect correction ...")
        # Plot PCA of uncorrected imputed data
        # transpose the dataframe to get samples as rows and features as columns
        self.df_pca = compute_pca(
            self.df.T,
            n_components=(n_components if n_components else round(len(set(self.batch_index)) / 3)),
        )

        # add batch information to the dataframe
        self.df_pca["batch"] = self.df_pca.index.str.split("-").str[0]

        # plot PC1 vs PC2 with batch information using seaborn
        # put the legend outside the plot
        # save the plot as a png file
        plot_pca(
            self.df_pca,
            title="PCA plot of uncorrected data",
            output_file="pca_uncorrected.png",
        )

        # keep samples only in tissue_part from metadata
        # TODO: specify covariates
        if tissue_parts_to_keep:
            self.metadata = self.metadata[self.metadata["tissue_part"].isin(tissue_parts_to_keep)]
            samples_to_keep = self.metadata["sample_id"].tolist()

            # keep samples in df that are also in samples_to_keep
            self.df = self.df[[s for s in self.df.columns if s in samples_to_keep]]

        # 2. Apply batch correction with covariate information
        # Before apply batch correction, filter out batches with just one sample (otherwise the batch correction will fail).
        batch_index = get_batch_info_from_sample_names(self.df.columns.tolist())
        self.df = remove_single_sample_batches(self.df, batch_index)

        # get covariate information from metadata.
        columns = self.df.columns.tolist()
        self.metadata = self.metadata[self.metadata["sample_id"].isin(columns)]
        # reorder metadata to match the order of columns in df
        self.metadata = self.metadata.reset_index(drop=True)
        self.metadata = self.metadata.set_index("sample_id").reindex(columns, axis=0).reset_index()
        if self.covariate:
            # get the covariates from metadata as a list
            covariates_index = self.metadata[self.covariate].tolist()
        else:
            covariates_index = []

        # apply batch correction
        self.df_corrected = apply_batch_correction(
            self.df, self.batch_index, covs=covariates_index
        )
        logger.info(self.df_corrected)

        # plot PCA of corrected data
        # transpose the dataframe to get samples as rows and features as columns
        # add batch information to the dataframe
        self.df_pca["batch"] = self.df_pca.index.str.split("-").str[0]

        # plot PC1 vs PC2 with batch information using seaborn
        # put the legend outside the plot
        # save the plot as a png file
        plot_pca(
            self.df_pca,
            title="PCA plot of corrected data",
            output_file="pca_corrected.png",
        )
