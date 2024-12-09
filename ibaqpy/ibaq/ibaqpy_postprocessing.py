import logging
import pandas as pd
from ibaqpy.ibaq.ibaqpy_commons import IBAQ, IBAQ_NORMALIZED, IBAQ_PPB, IBAQ_LOG, TPA, COPYNUMBER


def remove_samples_low_protein_number(ibaq_df: pd.DataFrame, min_protein_num: int) -> pd.DataFrame:
    """
    This functions takes an ibaq Dataframe with the following columns:
    - ProteinName
    - SampleID
    - Condition
    - iBAQ
    ...
    ...
    ...
    - iBAQLog

    and removes samples with low protein number. The minimum number of proteins in a sample is defined by the min_protein_num parameter.

    :param ibaq_df: pd.DataFrame
    :param min_protein_num: int
    :return: pd.DataFrame
    """

    # Count the number of proteins per sample
    protein_num = ibaq_df.groupby("SampleID")["ProteinName"].nunique()

    # Get the samples with more than min_protein_num proteins
    samples_to_keep = protein_num[protein_num >= min_protein_num].index
    samples_to_remove = protein_num[protein_num < min_protein_num].index

    logging.info(
        "The number of samples with number of proteins lower than {} is {}".format(
            min_protein_num, len(samples_to_remove)
        )
    )

    # Filter the samples
    ibaq_df = ibaq_df[ibaq_df["SampleID"].isin(samples_to_keep)]
    return ibaq_df


def remove_missing_values(
    ibaq_df: pd.DataFrame, missingness_percentage: float = 0.3, expression_column: str = IBAQ
) -> pd.DataFrame:
    """
    This functions takes an ibaq Dataframe with the following columns:
    - ProteinName
    - SampleID
    - Condition
    - iBAQ
    ...
    ...
    ...
    - iBAQLog

    and removes the samples that for all the proteins has missing values for a given expression column
    and remove the samples that have a percentage of missing values higher than missingness_percentage.

    :param ibaq_df: pd.DataFrame
    :param missingness_percentage: float
    :param expression_column: str The expression column to consider missigness percentage
    :return: pd.DataFrame
    """

    # Get the missingness percentage per sample
    missingness = ibaq_df.groupby("SampleID")[expression_column].apply(
        lambda x: x.isnull().sum() / len(x)
    )

    # Get the samples with missingness percentage lower than missingness_percentage
    samples_to_keep = missingness[missingness <= missingness_percentage].index
    samples_to_remove = missingness[missingness > missingness_percentage].index

    logging.info(
        "The number of samples with missingness percentage higher than {} is {}".format(
            missingness_percentage, len(samples_to_remove)
        )
    )

    # Filter the samples
    ibaq_df = ibaq_df[ibaq_df["SampleID"].isin(samples_to_keep)]
    return ibaq_df


def describe_expression_metrics(
    ibaq_df: pd.DataFrame
) -> pd.DataFrame:
    """
    This function takes an ibaq Dataframe with the following columns:
    - ProteinName
    - SampleID
    - Condition
    - iBAQ
    ...
    ...
    ...
    - iBAQLog

    and return a dataframe for each expression column with the following metrics:
    - standard deviation across samples
    - mean across samples
    - minimum across samples
    - 25th percentile across samples
    - median across samples
    - 75th percentile across samples
    - maximum across samples

    :param ibaq_df: pd.DataFrame
    :return: pd.DataFrame
    """

    possible_expresssion_values = [IBAQ, IBAQ_NORMALIZED, IBAQ_LOG, IBAQ_PPB, TPA, COPYNUMBER]

    # Define the expression columns
    expression_columns = [col for col in ibaq_df.columns if col in possible_expresssion_values]

    # Get the metrics
    metrics = ibaq_df.groupby("ProteinName")[expression_columns].describe()
    return metrics


