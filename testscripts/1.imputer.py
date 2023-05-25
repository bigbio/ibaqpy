"""
Script to perform imputation of missing values using KNNImputer.

"""

# import libraries
from typing import List, Optional, Union
from sklearn.impute import KNNImputer
import pandas as pd


def impute_missing_values(
    data: Optional[Union[pd.DataFrame, List[pd.DataFrame]]],
    n_neighbors=5,
    weights="uniform",
    metric="nan_euclidean",
    keep_empty_features=True,
) -> Union[pd.DataFrame, List[pd.DataFrame]]:
    """
    Impute missing values in a DataFrame or each DataFrame in a list using KNNImputer.

    Parameters
    ----------
    data : Union[pd.DataFrame, List[pd.DataFrame]]
        A pandas DataFrame or list of pandas DataFrames with missing values.
    n_neighbors : int, optional
        Number of neighboring samples to use for imputation. Default is 5.
    weights : str, optional
        Weight function used in prediction. Default is 'uniform'.
    metric : str, optional
        Distance metric for searching neighbors. Default is 'nan_euclidean'.
    keep_empty_features : bool, optional
        Whether to keep empty features (no known samples). Default is True.

    Returns
    -------
    Union[pd.DataFrame, List[pd.DataFrame]]
        A pandas DataFrame or list of pandas DataFrames with imputed missing values.
    """
    imputer = KNNImputer(
        n_neighbors=n_neighbors,
        weights=weights,
        metric=metric,
        keep_empty_features=keep_empty_features,
    )

    if isinstance(data, pd.DataFrame):
        # If it's a single DataFrame, transform it and return immediately
        return pd.DataFrame(imputer.fit_transform(data),
                            columns=data.columns,
                            index=data.index)
    else:
        # Otherwise, use list comprehension to apply the imputer to each DataFrame
        return [
            pd.DataFrame(imputer.fit_transform(t), columns=t.columns, index=t.index)
            for t in data
        ]


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
        raise ValueError(
            f"'{cov_index_col}' is not a column in the provided DataFrame."
        )

    # Use list comprehension to create the list of dataframes
    df_split = [df_group for _, df_group in df.groupby(cov_index_col)]

    return df_split


def filter_missing_value_by_group(df_input, cov_index_col, non_missing_percent_to_keep):
    """
    Filters the dataframe by keeping columns with at least a specified percent of non-missing values
    in each covariate_index group.

    Parameters:
    df_input (pandas.DataFrame): The input dataframe.
    cov_index_col (str): The name of the column to group by.
    non_missing_percent_to_keep (float): The minimum percentage of non-missing values to keep a column.

    Returns:
    pandas.DataFrame: The filtered dataframe.
    """
    return df_input.groupby(cov_index_col).apply(
        lambda x: x.dropna(thresh=non_missing_percent_to_keep * len(x), axis=1)
    )


# import data with first column as index
path = "/Users/enrique/projects/local/heart-proteomics"
df = pd.read_table(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm.tsv", index_col=0
)

# import metadata file for sample annotation (covariates index)
df_metadata = pd.read_csv(
    f"{path}/data/meta_data/heart_samples_meta_data_with_covs.tsv",
    sep="\t",
    index_col=False,
)

# transpose dataframe to have samples as rows and proteins as columns (a.k.a. features).
# keep index as column names
df = df.T

# Keep only columns 'sample_id' and 'covariate_index' from df_metadata
df_metadata = df_metadata[["sample_id", "covariate_index"]]

# Set index with column 'sample_id'
df_metadata.set_index("sample_id", inplace=True)

# merge df and df_metadata (annotate samples with covariate_index)
df = df.join(df_metadata, how="left")

# Keep only rows with covariate_index within 12, 11, 4, 6.
groups_to_keep = [12, 11, 4, 6]
df = df[df["covariate_index"].isin(groups_to_keep)]

print(df.head())

# keep columns with at least 30% of non-missing values in each covariate_index group
df = filter_missing_value_by_group(df,
                                   cov_index_col="covariate_index",
                                   non_missing_percent_to_keep=0.3)

# split df by covariate_index
df_list = split_df_by_column(df, cov_index_col="covariate_index")

# impute missing values with KNNImputer for every df in df_list
df_list = impute_missing_values(df_list)

# concatenate all dataframes in df_list into one dataframe
df = pd.concat(df_list)

# drop column 'covariate_index'
df.drop("covariate_index", axis=1, inplace=True)

# transpose dataframe to have proteins as rows and samples as columns as original.
# keep column names from index
df = df.T

# save dataframe as a tsv file
df.to_csv(
    f"{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed.tsv",
    sep="\t",
)
