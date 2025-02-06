from typing import Optional, Union, List

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer, SimpleImputer

from ibaqpy.ibaq.ibaqpy_commons import PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION, NORM_INTENSITY


def impute_missing_values(
    data: Optional[Union[pd.DataFrame, List[pd.DataFrame], None]],
    method: str = "knn",
    n_neighbors: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
    keep_empty_features: bool = True,
    fill_value: float = 0.0,
) -> Union[pd.DataFrame, List[pd.DataFrame], None]:
    """
    Impute missing values in a DataFrame or a list of DataFrames using KNN, mean, median, most frequent, or a specific value.

    Parameters
    ----------
    data : Optional[Union[pd.DataFrame, List[pd.DataFrame]]]
        A pandas DataFrame or a list of pandas DataFrames containing missing values to be imputed.
        The DataFrame(s) must adhere to the following format:
        - Rows represent samples (observations).
        - Columns represent features (variables).
        - Contain only numerical columns (e.g., float or int).
        - Missing values must be explicitly represented as `np.nan` or `pd.NA`.
        - Columns with non-numerical data (e.g., categorical or text) should be preprocessed
          (e.g., encoded into numerical values) before using this function.
        - Features (columns) with entirely missing values are handled based on the
          `keep_empty_features` parameter.
    method : str, optional
        The imputation method to use. Options are:
        - "knn" (default): Use K-Nearest Neighbors imputation.
        - "mean": Impute using the mean of each column.
        - "median": Impute using the median of each column.
        - "most_frequent": Impute using the most frequent value of each column.
        - "constant": Impute using a specific value provided via `fill_value`.
    n_neighbors : int, optional
        The number of neighboring samples to use for KNN imputation. Default is 5.
    weights : str, optional
        The weight function used in KNN prediction. Can be 'uniform' or 'distance'. Default is 'uniform'.
    metric : str, optional
        The distance metric used for finding neighbors in KNN. Default is 'nan_euclidean'.
    fill_value : float, optional
        The constant value to use for imputation when `method` is "constant". Default is 0.0.
    keep_empty_features : bool, optional
        Whether to keep features that are entirely empty (i.e., all values are NaN). Default is True.

    Returns
    -------
    Union[pd.DataFrame, List[pd.DataFrame]]
        A pandas DataFrame or a list of pandas DataFrames with imputed missing values.
        If the input is None, the function will return None.

    Notes
    -----
    - This function uses sklearn's KNNImputer and SimpleImputer for imputing missing values.
    - The `nan_euclidean` metric is specifically designed to handle NaN values during distance computation.
    - Column names and indices are preserved in the output.
    - Ensure the input data is numerical and properly formatted for the imputer.
    """
    if data is None:
        # placeholder for further implementation
        return None

    if method not in {"knn", "mean", "median", "constant", "most_frequent"}:
        raise ValueError(
            "Invalid method. Choose from 'knn', 'mean', 'median', 'most_frequent', or 'constant'."
        )

    if method == "knn":
        imputer = KNNImputer(
            n_neighbors=n_neighbors,
            weights=weights,
            metric=metric,
            keep_empty_features=keep_empty_features,
        )
    else:
        strategy = method
        imputer = SimpleImputer(strategy=strategy, fill_value=fill_value)

    def impute(df: pd.DataFrame) -> pd.DataFrame:
        imputed_data = imputer.fit_transform(df)
        return pd.DataFrame(imputed_data, columns=df.columns, index=df.index)

    if isinstance(data, pd.DataFrame):
        # Impute missing values for a single DataFrame
        return impute(data)
    elif isinstance(data, list) and all(isinstance(df, pd.DataFrame) for df in data):
        # Impute missing values for a list of DataFrames
        return [impute(df) for df in data]
    else:
        raise ValueError(
            "The input data must be a pandas DataFrame, a list of DataFrames, or None."
        )
