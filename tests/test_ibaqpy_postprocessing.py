from pathlib import Path
import pandas as pd

from ibaqpy.ibaq.ibaqpy_commons import SAMPLE_ID, IBAQ_NORMALIZED
from ibaqpy.ibaq.ibaqpy_postprocessing import (
    remove_samples_low_protein_number,
    remove_missing_values,
    describe_expression_metrics,
)
import logging

TESTS_DIR = Path(__file__).parent


def test_remove_samples_low_protein_number():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    ibaq_test = TESTS_DIR / "example/PXD017834-example-ibaq.tsv"
    ibaq_df = pd.read_csv(ibaq_test, sep="\t")
    number_samples = len(ibaq_df[SAMPLE_ID].unique())
    logging.info("The number of samples in the dataframe {}".format(number_samples))

    new_ibaq = remove_samples_low_protein_number(ibaq_df, min_protein_num=286)

    number_samples = len(new_ibaq[SAMPLE_ID].unique())
    logging.info(
        "The number of samples with number of proteins higher than 286 is {}".format(
            number_samples
        )
    )


def test_remove_missing_values():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    ibaq_test = TESTS_DIR / "example/PXD017834-example-ibaq.tsv"
    ibaq_df = pd.read_csv(ibaq_test, sep="\t")
    number_samples = len(ibaq_df[SAMPLE_ID].unique())
    logging.info("The number of samples in the dataframe {}".format(number_samples))
    new_ibaq = remove_missing_values(
        ibaq_df, missingness_percentage=1, expression_column=IBAQ_NORMALIZED
    )
    number_samples = len(new_ibaq[SAMPLE_ID].unique())
    logging.info(
        "The number of samples with less than 1% of missing values is {}".format(number_samples)
    )


def test_describe_expression_metrics():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    ibaq_test = TESTS_DIR / "example/PXD017834-example-ibaq.tsv"
    ibaq_df = pd.read_csv(ibaq_test, sep="\t")

    metrics = describe_expression_metrics(ibaq_df)
    logging.info(metrics)
