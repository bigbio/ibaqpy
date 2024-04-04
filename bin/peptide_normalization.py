#!/usr/bin/env python
import gc
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import os
import random
import uuid
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame
import pyarrow.parquet as pq
from parquet import Feature
from normalize_methods import normalize,normalize_run

from ibaq.ibaqpy_commons import (
    BIOREPLICATE,
    CHANNEL,
    CONDITION,
    PARQUET_COLUMNS,
    FRACTION,
    FRAGMENT_ION,
    INTENSITY,
    ISOTOPE_LABEL_TYPE,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PEPTIDE_CHARGE,
    PEPTIDE_SEQUENCE,
    PROTEIN_NAME,
    REFERENCE,
    RUN,
    SAMPLE_ID,
    STUDY_ID,
    TMT16plex,
    TMT11plex,
    TMT10plex,
    TMT6plex,
    ITRAQ4plex,
    ITRAQ8plex,
    get_canonical_peptide,
    get_spectrum_prefix,
    get_study_accession,
    parquet_map,
    parse_uniprot_accession,
    plot_box_plot,
    plot_distributions,
    remove_contaminants_entrapments_decoys,
    remove_protein_by_ids,
    sum_peptidoform_intensities,
    get_peptidoform_normalize_intensities,
    average_peptide_intensities,
    print_help_msg,
)


def print_dataset_size(dataset: DataFrame, message: str, verbose: bool) -> None:
    if verbose:
        print(message + str(len(dataset.index)))

def recover_df(df):
    """
    This function is aimed to recover data shape.
    """
    samples = df.columns.tolist()
    out = pd.DataFrame()
    for sample in samples:
        samples_df = df[sample].dropna()
        samples_df = samples_df.reset_index()
        samples_df['SampleID'] = sample
        samples_df.rename(columns={
            sample:NORM_INTENSITY
        },inplace=True)
        out = pd.concat([out,samples_df])
    out.reset_index(inplace=True,drop=True)
    return out

def analyse_sdrf(sdrf_path: str, compression: bool) -> tuple:
    """
    This function is aimed to parse SDRF and return four objects:
    1. sdrf_df: A dataframe with channels and references annoted.
    2. label: Label type of the experiment. LFQ, TMT or iTRAQ.
    3. sample_names: A list contains all sample names.
    4. choice: A dictionary caontains key-values between channel
        names and numbers.
    :param sdrf_path: File path of SDRF.
    :param compression: Whether compressed.
    :return:
    """
    sdrf_df = pd.read_csv(sdrf_path, sep="\t", compression=compression)
    sdrf_df[REFERENCE] = sdrf_df["comment[data file]"].apply(get_spectrum_prefix)

    labels = set(sdrf_df["comment[label]"])
    # Determine label type
    label, choice = get_label(labels)
    if label == "TMT":
        choice_df = (
            pd.DataFrame.from_dict(choice, orient="index", columns=[CHANNEL])
            .reset_index()
            .rename(columns={"index": "comment[label]"})
        )
        sdrf_df = sdrf_df.merge(choice_df, on="comment[label]", how="left")
    elif label == "ITRAQ":
        choice_df = (
            pd.DataFrame.from_dict(choice, orient="index", columns=[CHANNEL])
            .reset_index()
            .rename(columns={"index": "comment[label]"})
        )
        sdrf_df = sdrf_df.merge(choice_df, on="comment[label]", how="left")
    sample_names = sdrf_df["source name"].unique().tolist()

    return sdrf_df, label, sample_names, choice


def analyse_feature_df(feature_df: pd.DataFrame) -> tuple:
    """Return label type, sample names and choice dict by iterating parquet.

    :param parquet_path: Feature parquet path.
    :param batch_size: Iterate batch size, defaults to 100000
    :return: Label type, sample names and choice dict
    """
    samples = feature_df["sample_accession"].unique().tolist()
    labels = feature_df["isotope_label_type"].unique().tolist()
    # Determine label type
    label, choice = get_label(labels)

    return label, samples, choice


def analyse_feature_parquet(parquet_path: str, batch_size: int = 100000) -> tuple:
    """Return label type, sample names and choice dict by iterating parquet.

    :param parquet_path: Feature parquet path.
    :param batch_size: Iterate batch size, defaults to 100000
    :return: Label type, sample names and choice dict
    """
    parquet_chunks = read_large_parquet(parquet_path, batch_size)
    labels, samples = list(), list()
    for chunk in parquet_chunks:
        samples.extend(chunk["sample_accession"].unique().tolist())
        labels.extend(chunk["isotope_label_type"].unique().tolist())
        samples = list(set(samples))
        labels = list(set(labels))
    # Determine label type
    label, choice = get_label(labels)

    return label, samples, choice


def read_large_parquet(parquet_path: str, batch_size: int = 100000):
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batch_df = batch.to_pandas()
        yield batch_df


def get_label(labels: list) -> (str, dict):
    """Return label type and choice dict according to labels list.

    :param labels: Labels from SDRF.
    :return: Tuple contains label type and choice dict.
    """
    choice = None
    if len(labels) == 1:
        label = "LFQ"
    elif "TMT" in ",".join(labels) or "tmt" in ",".join(labels):
        if (
            len(labels) > 11
            or "TMT134N" in labels
            or "TMT133C" in labels
            or "TMT133N" in labels
            or "TMT132C" in labels
            or "TMT132N" in labels
        ):
            choice = TMT16plex
        elif len(labels) == 11 or "TMT131C" in labels:
            choice = TMT11plex
        elif len(labels) > 6:
            choice = TMT10plex
        else:
            choice = TMT6plex
        label = "TMT"
    elif "ITRAQ" in ",".join(labels) or "itraq" in ",".join(labels):
        if len(labels) > 4:
            choice = ITRAQ8plex
        else:
            choice = ITRAQ4plex
        label = "ITRAQ"
    else:
        exit("Warning: Only support label free, TMT and ITRAQ experiment!")
    return label, choice


def msstats_common_process(data_df: pd.DataFrame) -> pd.DataFrame:
    """Apply common process on data.

    :param data_df: Feature data in dataframe.
    :return: Processed data.
    """
    data_df.rename(
        columns={
            "ProteinName": PROTEIN_NAME,
            "PeptideSequence": PEPTIDE_SEQUENCE,
            "PrecursorCharge": PEPTIDE_CHARGE,
            "Run": RUN,
            "Condition": CONDITION,
            "Intensity": INTENSITY,
        },
        inplace=True,
    )
    data_df[REFERENCE] = data_df[REFERENCE].apply(get_spectrum_prefix)

    return data_df


def parquet_common_process(
    data_df: pd.DataFrame, label: str, choice: dict
) -> pd.DataFrame:
    """Apply common process on data.

    :param data_df: Feature data in dataframe.
    :return: Processed data.
    """
    data_df = data_df.rename(columns=parquet_map)
    data_df[PROTEIN_NAME] = data_df.apply(lambda x: ",".join(x[PROTEIN_NAME]), axis=1)
    if label == "LFQ":
        data_df.drop(CHANNEL, inplace=True, axis=1)
    else:
        data_df[CHANNEL] = data_df[CHANNEL].map(choice)

    return data_df


def merge_sdrf(
    label: str, sdrf_df: pd.DataFrame, data_df: pd.DataFrame
) -> pd.DataFrame:
    if label == "LFQ":
        result_df = pd.merge(
            data_df,
            sdrf_df[["source name", REFERENCE]],
            how="left",
            on=[REFERENCE],
        )
    elif label == "TMT":
        result_df = pd.merge(
            data_df,
            sdrf_df[["source name", REFERENCE, CHANNEL]],
            how="left",
            on=[REFERENCE, CHANNEL],
        )
    elif label == "ITRAQ":
        result_df = pd.merge(
            data_df,
            sdrf_df[["source name", REFERENCE, CHANNEL]],
            how="left",
            on=[REFERENCE, CHANNEL],
        )
    result_df.rename(columns={"source name": SAMPLE_ID}, inplace=True)
    result_df = result_df[result_df["Condition"] != "Empty"]

    return result_df


def data_common_process(data_df: pd.DataFrame, min_aa: int) -> pd.DataFrame:
    # Remove 0 intensity signals from the data
    data_df = data_df[data_df[INTENSITY] > 0]
    data_df = data_df[data_df["Condition"] != "Empty"]

    def map_canonical_seq(data_df: pd.DataFrame) -> (pd.DataFrame, dict):
        modified_seqs = data_df[PEPTIDE_SEQUENCE].unique().tolist()
        canonical_seqs = [get_canonical_peptide(i) for i in modified_seqs]
        inner_canonical_dict = dict(zip(modified_seqs, canonical_seqs))
        data_df[PEPTIDE_CANONICAL] = data_df[PEPTIDE_SEQUENCE].map(inner_canonical_dict)

        return data_df, inner_canonical_dict

    if PEPTIDE_CANONICAL not in data_df.columns:
        data_df, inner_canonical_dict = map_canonical_seq(data_df)
        data_df[PEPTIDE_CANONICAL] = data_df[PEPTIDE_SEQUENCE].map(inner_canonical_dict)
    # Filter peptides with less amino acids than min_aa (default: 7)
    data_df = data_df[
        data_df.apply(lambda x: len(x[PEPTIDE_CANONICAL]) >= min_aa, axis=1)
    ]
    data_df[PROTEIN_NAME] = data_df[PROTEIN_NAME].apply(parse_uniprot_accession)
    #data_df[STUDY_ID] = data_df[SAMPLE_ID].apply(get_study_accession)
    if FRACTION not in data_df.columns:
        data_df[FRACTION] = 1
        data_df = data_df[
            [
                PROTEIN_NAME,
                PEPTIDE_SEQUENCE,
                PEPTIDE_CANONICAL,
                PEPTIDE_CHARGE,
                INTENSITY,
                CONDITION,
                RUN,
                BIOREPLICATE,
                FRACTION,
                ISOTOPE_LABEL_TYPE,
                SAMPLE_ID,
            ]
        ]
    data_df[CONDITION] = pd.Categorical(data_df[CONDITION])
    #data_df[STUDY_ID] = pd.Categorical(data_df[STUDY_ID])
    data_df[SAMPLE_ID] = pd.Categorical(data_df[SAMPLE_ID])

    return data_df

def intensity_normalization(
    dataset: DataFrame,
    field: str,
    class_field: str,
    scaling_method: str = "quantile",
) -> DataFrame:
    cols_to_keep = [
        PROTEIN_NAME,
        PEPTIDE_CANONICAL,
        PEPTIDE_SEQUENCE,
        PEPTIDE_CHARGE,
        SAMPLE_ID,
        BIOREPLICATE,
        CONDITION,
        NORM_INTENSITY,
    ]
    # TODO add imputation and/or removal to those two norm strategies
    if scaling_method == "msstats":
        # For TMT normalization
        if "Channel" in dataset.columns:
            g = dataset.groupby(["Run", "Channel"])[field].apply(np.nanmedian)
            g.name = "RunMedian"
            dataset = dataset.join(g, on=["Run", "Channel"])
            median_baseline = dataset.drop_duplicates(subset=["Run", "Channel", field])[
                field
            ].median()
            dataset[NORM_INTENSITY] = (
                dataset[field] - dataset["RunMedian"] + median_baseline
            )
        else:
            g = dataset.groupby(["Run", "Fraction"])[field].apply(np.nanmedian)
            g.name = "RunMedian"
            dataset = dataset.join(g, on=["Run", "Fraction"])
            dataset["FractionMedian"] = (
                dataset["RunMedian"].groupby(dataset["Fraction"]).transform("median")
            )
            dataset[NORM_INTENSITY] = (
                dataset[field] - dataset["RunMedian"] + dataset["FractionMedian"]
            )
        return dataset[cols_to_keep]

    else:
        # pivot to have one col per sample
        print("Transforming to wide format dataset size {}".format(len(dataset.index)))
        normalize_df = pd.pivot_table(
            dataset,
            index=[
                PEPTIDE_SEQUENCE,
                PEPTIDE_CANONICAL,
                PEPTIDE_CHARGE,
                FRACTION,
                RUN,
                BIOREPLICATE,
                PROTEIN_NAME,
                STUDY_ID,
                CONDITION,
            ],
            columns=class_field,
            values=field,
            aggfunc={field: np.nanmean},
            observed=True,
        )
        normalize_df = normalize(normalize_df,scaling_method)
        # TODO: When restoring the pivot table here, the previous grouping caused
        # the dataframe to produce a large number of rows with NORM_INTENSITY of
        # NA at melt. This results in an unbearable memory consumption.

        normalize_df = recover_df(normalize_df)
        normalize_df = normalize_df.drop_duplicates()
        print(normalize_df.head())
        return normalize_df[cols_to_keep]


def remove_low_frequency_peptides_(
    dataset_df: DataFrame, percentage_samples: float = 0.20
):
    """
    Remove peptides that are present in less than 20% of the samples.
    :param dataset_df: dataframe with the data
    :param percentage_samples: percentage of samples
    :return:
    """
    #= dataset_df[[SAMPLE_ID, CONDITION]].drop_duplicates(subset=[SAMPLE_ID])
    c_map = dataset_df[[SAMPLE_ID, CONDITION]].drop_duplicates(subset=[SAMPLE_ID]).set_index(SAMPLE_ID).to_dict()[CONDITION]
    normalize_df = pd.pivot_table(
        dataset_df,
        index=[PEPTIDE_CANONICAL, PROTEIN_NAME],
        columns=SAMPLE_ID,
        values=NORM_INTENSITY,
        aggfunc={NORM_INTENSITY: np.nanmean},
        observed=True,
    )
    del dataset_df
    gc.collect()
    # Count the number of null values in each row
    null_count = normalize_df.isnull().sum(axis=1)
    # Find the rows that have null values above the threshold
    rows_to_drop = null_count[
        null_count >= (1 - percentage_samples) * normalize_df.shape[1]
    ].index
    # Drop the rows with too many null values
    normalize_df = normalize_df.drop(rows_to_drop)

    # Remove rows with non-null values in only one column
    normalize_df = normalize_df[
        normalize_df.notnull().sum(axis=1) != 1
    ]
    normalize_df = recover_df(normalize_df)
    # recover condition column
    normalize_df.loc[:,CONDITION] =  normalize_df[SAMPLE_ID].map(c_map)

    # Remove rows with null values in NORMALIZE_INTENSITY
    normalize_df.dropna(subset=[NORM_INTENSITY], inplace=True)

    print(normalize_df.head())
    return normalize_df


def peptide_intensity_normalization(
    dataset_df: DataFrame, field: str, class_field: str, scaling_method: str
):
    """
    Normalize the peptide intensities using different methods.
    :param dataset_df: dataframe with the data
    :param field: field to normalize
    :param class_field: field to use as class
    :param scaling_method: method to use for the normalization
    :return:
    """
    # pivot to have one col per sample
    normalize_df = pd.pivot_table(
        dataset_df,
        index=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION],
        columns=class_field,
        values=field,
        aggfunc={field: np.nanmean},
        observed=True,
    )
    # need nomalize?
    normalize_df = recover_df(normalize_df)
    return normalize_df

@click.command()
@click.option(
    "-m", "--msstats", help="MsStats file import generated by quantms", default=None
)
@click.option(
    "-p", "--parquet", help="Parquet file import generated by quantms.io", default=None
)
@click.option(
    "-s", "--sdrf", help="SDRF file import generated by quantms", default=None
)
@click.option("--stream", help="Stream processing normalization", is_flag=True)
@click.option(
    "--chunksize",
    help="The number of rows of MSstats or parquet read using pandas streaming",
    default=1000000,
)
@click.option(
    "--min_aa", help="Minimum number of amino acids to filter peptides", default=7
)
@click.option(
    "--min_unique",
    help="Minimum number of unique peptides to filter proteins",
    default=2,
)
@click.option(
    "--remove_ids",
    help="Remove specific protein ids from the analysis using a file with one id per line",
)
@click.option(
    "--remove_decoy_contaminants",
    help="Remove decoy and contaminants proteins from the analysis",
    is_flag=True,
    default=False,
)
@click.option(
    "--remove_low_frequency_peptides",
    help="Remove peptides that are present in less than 20% of the samples",
    is_flag=True,
    default=False,
)
@click.option(
    "--output",
    help="Peptide intensity file including other all properties for normalization",
)
@click.option(
    "--skip_normalization", help="Skip normalization step", is_flag=True, default=False
)
@click.option(
    "--nmethod",
    help="Normalization method used to normalize intensities for all samples (options: quantile, msstats, qnorm)",
    default="quantile",
)
@click.option(
    "--pnormalization",
    help="Normalize the peptide intensities using different methods (options: qnorm)",
    is_flag=True,
)
@click.option(
    "--compress",
    help="Read the input peptides file in compress gzip file",
    is_flag=True,
)
@click.option(
    "--violin",
    help="Use violin plot instead of boxplot for distribution representations",
    is_flag=True,
)
@click.option(
    "--verbose",
    help="Print addition information about the distributions of the intensities, number of peptides remove "
    "after normalization, etc.",
    is_flag=True,
)
@click.option(
    "--qc_report",
    help="PDF file to store multiple QC images",
    default="peptideNorm-QCprofile.pdf",
)
def peptide_normalization(
    msstats: str,
    parquet: str,
    sdrf: str,
    stream: bool,
    chunksize: int,
    min_aa: int,
    min_unique: int,
    remove_ids: str,
    remove_decoy_contaminants: bool,
    remove_low_frequency_peptides: bool,
    output: str,
    skip_normalization: bool,
    nmethod: str,
    pnormalization: bool,
    compress: bool,
    violin: bool,
    verbose: bool,
    qc_report: str,
) -> None:
    if output is None:
        print_help_msg(peptide_normalization)
        exit(1)

    if parquet is None and (msstats is None or sdrf is None):
        print_help_msg(peptide_normalization)
        exit(1)

    if pnormalization and nmethod not in ["qnorm", "quantile"]:
        exit(
            "Peptide intensity normalization works only with qnorm or quantile methods!"
        )

    if verbose:
        log_after_norm = not log2

    pd.set_option("display.max_columns", None)
    compression_method = "gzip" if compress else None
    print("Loading data..")
    '''
    if parquet is None:
        # Read the msstats file
        feature_df = pd.read_csv(
            msstats,
            sep=",",
            compression=compression_method,
            dtype={CONDITION: "category", ISOTOPE_LABEL_TYPE: "category"},
        )

        # Read the sdrf file
        sdrf_df, label, sample_names, choice = analyse_sdrf(
            sdrf, compression_method
        )
        print(sdrf_df)

        # Merged the SDRF with the Resulted file
        dataset_df = msstats_common_process(feature_df)
        dataset_df = merge_sdrf(label, sdrf_df, feature_df)
        # Remove the intermediate variables and free the memory
        del feature_df, sdrf_df
        gc.collect()
    else:
    '''
    F = Feature(parquet)
    #dataset_df = pd.read_parquet(parquet,columns=PARQUET_COLUMNS)
    header= False
    for samples,df in F.iter_samples():
        for sample in samples:
            dataset_df = df[df['sample_accession']==sample].copy()
            dataset_df = dataset_df[PARQUET_COLUMNS]
            label, sample_names, choice = analyse_feature_df(dataset_df)
            dataset_df = parquet_common_process(dataset_df, label, choice)
            dataset_df = data_common_process(dataset_df, min_aa)
            # Only proteins with unique peptides number greater than min_unique (default: 2) are retained
            unique_peptides = set(
                dataset_df.groupby(PEPTIDE_CANONICAL)
                .filter(lambda x: len(set(x[PROTEIN_NAME])) == 1)[PEPTIDE_CANONICAL]
                .tolist()
            )
            strong_proteins = set(
                dataset_df[dataset_df[PEPTIDE_CANONICAL].isin(unique_peptides)]
                .groupby(PROTEIN_NAME)
                .filter(lambda x: len(set(x[PEPTIDE_CANONICAL])) >= min_unique)[
                    PROTEIN_NAME
                ]
                .tolist()
            )
            dataset_df = dataset_df[dataset_df[PROTEIN_NAME].isin(strong_proteins)]

            print(f"Number of unique peptides: {len(unique_peptides)}")
            print(f"Number of strong proteins: {len(strong_proteins)}")

            print("Logarithmic if specified..")
            dataset_df.rename(columns={INTENSITY: NORM_INTENSITY},inplace=True)
            """     
            if log2:
                dataset_df[NORM_INTENSITY] = np.log2(dataset_df[NORM_INTENSITY]) 
            """

            # Remove high abundant and contaminants proteins and the outliers
            if remove_ids is not None:
                print("Remove proteins from file...")
                dataset_df = remove_protein_by_ids(dataset_df, remove_ids)
            if remove_decoy_contaminants:
                print("Remove decoy and contaminants...")
                dataset_df = remove_contaminants_entrapments_decoys(dataset_df)

            print_dataset_size(dataset_df, "Peptides after contaminants removal: ", verbose)
            print("Normalize intensities.. ")
            dataset_df.loc[:,NORM_INTENSITY] = np.log(dataset_df[NORM_INTENSITY])
            dataset_df = normalize_run(dataset_df,sdrf,nmethod)
            dataset_df.loc[:,NORM_INTENSITY] = normalize(dataset_df[NORM_INTENSITY],'max_min')
            print("Number of peptides after normalization: " + str(len(dataset_df.index)))
            print("Select the best peptidoform across fractions...")
            dataset_df = get_peptidoform_normalize_intensities(dataset_df)
            print(
                "Number of peptides after peptidofrom selection: "
                + str(len(dataset_df.index))
            )

            print("Sum all peptidoforms per Sample...")
            dataset_df = sum_peptidoform_intensities(dataset_df)
            print("Number of peptides after selection: " + str(len(dataset_df.index)))

            print("Average all peptidoforms per Peptide/Sample...")
            dataset_df = average_peptide_intensities(dataset_df)
            print("Number of peptides after average: " + str(len(dataset_df.index)))

            if remove_low_frequency_peptides and len(sample_names) > 1:
                print(dataset_df)
                dataset_df = remove_low_frequency_peptides_(dataset_df, 0.20)
                print_dataset_size(
                    dataset_df, "Peptides after remove low frequency peptides: ", verbose
                )


            print("Save the normalized peptide intensities...")
            if header:
                dataset_df.to_csv(output, index=False,header=False,mode='a+')
            else:
                dataset_df.to_csv(output, index=False)
                header = True
            


if __name__ == "__main__":
    peptide_normalization()
