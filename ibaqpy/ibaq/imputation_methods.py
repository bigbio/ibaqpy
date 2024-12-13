from typing import Optional, Union, List

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from ibaqpy.ibaq.ibaqpy_commons import PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION, NORM_INTENSITY


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
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)
    else:
        # Otherwise, use list comprehension to apply the imputer to each DataFrame
        return [
            pd.DataFrame(imputer.fit_transform(t), columns=t.columns, index=t.index) for t in data
        ]


def impute_missing_values(dataset_df, field, class_field):
    """
    Impute the missing values using different methods.
    :param dataset_df: dataframe with the data
    :param field: field to impute
    :param class_field: field to use as class
    :return:
    """
    normalize_df = pd.DataFrame()
    # group by condition to detect missing values
    for c, g in dataset_df.groupby(CONDITION):
        # pivot to have one col per sample
        group_normalize_df = pd.pivot_table(
            g,
            index=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION],
            columns=class_field,
            values=field,
            aggfunc={field: np.nanmean},
            observed=True,
        )

        # no missing values group -> only one sample
        if len(group_normalize_df.columns) < 2:
            group_normalize_df = group_normalize_df.reset_index()
            group_normalize_df = group_normalize_df.melt(
                id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION]
            )
            group_normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
            normalize_df = pd.concat([normalize_df, group_normalize_df], ignore_index=True)

    return normalize_df
