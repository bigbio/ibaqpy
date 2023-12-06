#!/usr/bin/env python
import gc
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qnorm
from scipy.stats import rankdata
import os
import random
import uuid
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame
import pyarrow.parquet as pq

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
    SEARCH_ENGINE,
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
    remove_contaminants_decoys,
    remove_protein_by_ids,
    sum_peptidoform_intensities,
    get_peptidoform_normalize_intensities,
    average_peptide_intensities,
    print_help_msg,
)


def print_dataset_size(dataset: DataFrame, message: str, verbose: bool) -> None:
    if verbose:
        print(message + str(len(dataset.index)))


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
    data_df[STUDY_ID] = data_df[SAMPLE_ID].apply(get_study_accession)
    if FRACTION not in data_df.columns:
        data_df[FRACTION] = 1
        data_df = data_df[
            [
                PROTEIN_NAME,
                PEPTIDE_SEQUENCE,
                PEPTIDE_CANONICAL,
                PEPTIDE_CHARGE,
                INTENSITY,
                REFERENCE,
                CONDITION,
                RUN,
                BIOREPLICATE,
                FRACTION,
                FRAGMENT_ION,
                ISOTOPE_LABEL_TYPE,
                STUDY_ID,
                SAMPLE_ID,
            ]
        ]
    data_df[CONDITION] = pd.Categorical(data_df[CONDITION])
    data_df[STUDY_ID] = pd.Categorical(data_df[STUDY_ID])
    data_df[SAMPLE_ID] = pd.Categorical(data_df[SAMPLE_ID])

    return data_df


# TODO: Here we present a method to apply quantile normalization. However,
# this function runs slowly and needs to be improved.
def quantile_normalize(
    data: pd.DataFrame, index: pd.core.indexes, cols: pd.core.indexes
) -> pd.DataFrame:
    """Quantile normalization which considers NA values.

    :param data: Pivot dataframe contains intensities for each peptides.
    :param index: Dataframe index.
    :param cols: Dataframe columns.
    :return: Dataframe normalizaed.
    """
    data = data.to_numpy()
    sorted_data = np.sort(data, axis=0)
    row_quantiles = np.nanmean(sorted_data, axis=1)
    rank_data = rankdata(data, axis=0, nan_policy="omit")
    for i in range(data.shape[1]):
        locs = (~np.isnan(data[:, i])).nonzero()
        data[locs, i] = row_quantiles[[int(i) for i in (rank_data[locs, i] - 1)[0]]]

    return pd.DataFrame(data, index=index, columns=cols)


def intensity_normalization(
    dataset: DataFrame,
    field: str,
    class_field: str,
    scaling_method: str = "quantile",
) -> DataFrame:
    # TODO add imputation and/or removal to those two norm strategies
    if scaling_method == "msstats":
        # For TMT normalization
        if "Channel" in dataset.columns:
            g = dataset.groupby(["Run", "Channel"])[field].apply(np.median)
            g.name = "RunMedian"
            dataset = dataset.join(g, on=["Run", "Channel"])
            median_baseline = dataset.drop_duplicates(subset=["Run", "Channel", field])[
                field
            ].median()
            dataset[NORM_INTENSITY] = (
                dataset[field] - dataset["RunMedian"] + median_baseline
            )
        else:
            g = dataset.groupby(["Run", "Fraction"])[field].apply(np.median)
            g.name = "RunMedian"
            dataset = dataset.join(g, on=["Run", "Fraction"])
            dataset["FractionMedian"] = (
                dataset["RunMedian"].groupby(dataset["Fraction"]).transform("median")
            )
            dataset[NORM_INTENSITY] = (
                dataset[field] - dataset["RunMedian"] + dataset["FractionMedian"]
            )
        return dataset

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
            aggfunc={field: np.mean},
            observed=True,
        )
        if scaling_method == "qnorm":
            normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        elif scaling_method == "quantile":
            normalize_df = quantile_normalize(
                normalize_df, normalize_df.index, normalize_df.columns
            )
        # TODO: When restoring the pivot table here, the previous grouping caused
        # the dataframe to produce a large number of rows with NORM_INTENSITY of
        # NA at melt. This results in an unbearable memory consumption.
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(
            id_vars=[
                PEPTIDE_SEQUENCE,
                PEPTIDE_CANONICAL,
                PEPTIDE_CHARGE,
                FRACTION,
                RUN,
                BIOREPLICATE,
                PROTEIN_NAME,
                STUDY_ID,
                CONDITION,
            ]
        )
        normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
        normalize_df = normalize_df[normalize_df[NORM_INTENSITY].notna()]
        normalize_df = normalize_df.drop_duplicates()
        print(normalize_df.head())
        return normalize_df

    return dataset


def remove_low_frequency_peptides_(
    dataset_df: DataFrame, percentage_samples: float = 0.20
):
    """
    Remove peptides that are present in less than 20% of the samples.
    :param dataset_df: dataframe with the data
    :param percentage_samples: percentage of samples
    :return:
    """

    normalize_df = pd.pivot_table(
        dataset_df,
        index=[PEPTIDE_CANONICAL, PROTEIN_NAME],
        columns=SAMPLE_ID,
        values=NORM_INTENSITY,
        aggfunc={NORM_INTENSITY: np.mean},
        observed=True,
    )
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
        normalize_df.notnull().sum(axis=1) != normalize_df.shape[1] - 1
    ]
    normalize_df = normalize_df.reset_index()
    normalize_df = normalize_df.melt(id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME])
    normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)

    # recover condition column
    normalize_df = normalize_df.merge(
        dataset_df[[SAMPLE_ID, CONDITION]].drop_duplicates(subset=[SAMPLE_ID]),
        on=SAMPLE_ID,
        how="left",
    )

    # Remove rows with null values in NORMALIZE_INTENSITY
    normalize_df = normalize_df[normalize_df[NORM_INTENSITY].notna()]

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
        aggfunc={field: np.mean},
        observed=True,
    )
    if scaling_method == "qnorm":
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
    elif scaling_method == "quantile":
        normalize_df = quantile_normalize(
            normalize_df, normalize_df.index, normalize_df.columns
        )
    else:
        exit(
            "Peptide intensity normalization works only with qnorm or quantile methods!"
        )
    normalize_df = normalize_df.reset_index()
    normalize_df = normalize_df.melt(
        id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION]
    )
    normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
    normalize_df = normalize_df[normalize_df[NORM_INTENSITY].notna()]
    return normalize_df


def impute_peptide_intensities(dataset_df, field, class_field):
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
            aggfunc={field: np.mean},
            observed=True,
        )

        # no missing values group -> only one sample
        if len(group_normalize_df.columns) < 2:
            group_normalize_df = group_normalize_df.reset_index()
            group_normalize_df = group_normalize_df.melt(
                id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION]
            )
            group_normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
            normalize_df = normalize_df.append(group_normalize_df, ignore_index=True)
        # else:
        #     # print ("nothing")
        #     # Impute the missing values
        #     # imputer = MissForest(max_iter=5)
        #     # imputed_data = imputer.fit_transform(group_normalize_df)
        #     # group_normalize_df = pd.DataFrame(imputed_data, columns=group_normalize_df.columns,
        #     #                                   index=group_normalize_df.index)
        #     # # Melt the dataframe
        #     # group_normalize_df = group_normalize_df.reset_index()
        #     # group_normalize_df = group_normalize_df.melt(id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION])
        #     # group_normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
        #     # normalize_df = normalize_df.append(group_normalize_df, ignore_index=True)

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
    "--log2",
    help="Transform to log2 the peptide intensity values before normalization",
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
    log2: bool,
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

    pd.set_option("display.max_columns", None)
    compression_method = "gzip" if compress else None
    print("Loading data..")

    if not stream:
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
            dataset_df = pd.read_parquet(parquet)[PARQUET_COLUMNS]
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
        dataset_df.loc[dataset_df.Intensity == 0, INTENSITY] = 1
        dataset_df[NORM_INTENSITY] = (
            np.log2(dataset_df[INTENSITY]) if log2 else dataset_df[INTENSITY]
        )
        dataset_df.drop(INTENSITY, axis=1, inplace=True)

        # Print the distribution of the original peptide intensities from quantms analysis
        if verbose:
            sample_names = set(dataset_df[SAMPLE_ID])
            plot_width = len(sample_names) * 0.5 + 10
            pdf = PdfPages(qc_report)
            density = plot_distributions(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=not log2,
                width=plot_width,
                title="Original peptidoform intensity distribution (no normalization)",
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=not log2,
                width=plot_width,
                title="Original peptidoform intensity distribution (no normalization)",
                violin=violin,
            )
            plt.show()
            pdf.savefig(box)

        # Remove high abundant and contaminants proteins and the outliers
        if remove_ids is not None:
            print("Remove proteins from file...")
            dataset_df = remove_protein_by_ids(dataset_df, remove_ids)
        if remove_decoy_contaminants:
            print("Remove decoy and contaminants...")
            dataset_df = remove_contaminants_decoys(dataset_df)

        print_dataset_size(dataset_df, "Peptides after contaminants removal: ", verbose)

        print("Normalize intensities.. ")
        # dataset_df = dataset_df.dropna(how="any")
        if not skip_normalization:
            dataset_df = intensity_normalization(
                dataset_df,
                field=NORM_INTENSITY,
                class_field=SAMPLE_ID,
                scaling_method=nmethod,
            )
        if verbose:
            log_after_norm = (
                nmethod == "msstats"
                or nmethod == "qnorm"
                or ((nmethod == "quantile" or nmethod == "robust") and not log2)
            )
            density = plot_distributions(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Peptidoform intensity distribution after normalization, method: "
                + nmethod,
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Peptidoform intensity distribution after normalization, method: "
                + nmethod,
                violin=violin,
            )
            plt.show()
            pdf.savefig(box)
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

        if verbose:
            log_after_norm = (
                nmethod == "msstats"
                or nmethod == "qnorm"
                or ((nmethod == "quantile" or nmethod == "robust") and not log2)
            )
            density = plot_distributions(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Peptide intensity distribution method: " + nmethod,
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Peptide intensity distribution method: " + nmethod,
                violin=violin,
            )
            plt.show()
            pdf.savefig(box)

        if remove_low_frequency_peptides and len(sample_names) > 1:
            print(dataset_df)
            dataset_df = remove_low_frequency_peptides_(dataset_df, 0.20)
            print_dataset_size(
                dataset_df, "Peptides after remove low frequency peptides: ", verbose
            )

        # Perform imputation using Random Forest in Peptide Intensities
        # TODO: Check if this is necessary (Probably we can do some research if imputation at peptide level is necessary
        # if impute:
        #     dataset_df = impute_peptide_intensities(dataset_df, field=NORM_INTENSITY, class_field=SAMPLE_ID)

        if pnormalization:
            print("Normalize at Peptide level...")
            dataset_df = peptide_intensity_normalization(
                dataset_df,
                field=NORM_INTENSITY,
                class_field=SAMPLE_ID,
                scaling_method=nmethod,
            )

        if verbose:
            log_after_norm = (
                nmethod == "msstats"
                or nmethod == "qnorm"
                or ((nmethod == "quantile" or nmethod == "robust") and not log2)
            )
            density = plot_distributions(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Normalization at peptide level method: " + nmethod,
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                dataset_df,
                NORM_INTENSITY,
                SAMPLE_ID,
                log2=log_after_norm,
                width=plot_width,
                title="Normalization at peptide level method: " + nmethod,
                violin=violin,
            )
            plt.show()
            pdf.savefig(box)
            pdf.close()

        print("Save the normalized peptide intensities...")
        dataset_df.to_csv(output, index=False, sep=",")
    else:
        if parquet is None:
            sdrf_df, label, sample_names, choice = analyse_sdrf(
                sdrf, compression_method
            )
            msstats_chunks = pd.read_csv(
                msstats,
                sep=",",
                compression=compression_method,
                dtype={CONDITION: "category", ISOTOPE_LABEL_TYPE: "category"},
                chunksize=chunksize,
            )
        else:
            label, sample_names, choice = analyse_feature_parquet(
                parquet, batch_size=chunksize
            )
            msstats_chunks = read_large_parquet(parquet, batch_size=chunksize)
        sample_number = len(sample_names)

        # TODO: Stream processing to obtain strong proteins with more than 2 uniqe peptides
        temp = f"Temp-{str(uuid.uuid4())}/"
        os.mkdir(temp)
        print(f"INFO: Writing files into {temp}...")
        unique_peptides = {}
        group_intensities = {}
        quantile = {}
        print("INFO: First iteration to get unique peptides and strong proteins...")
        for msstats_df in msstats_chunks:
            if parquet is None:
                msstats_df = msstats_common_process(msstats_df)
                msstats_df = merge_sdrf(label, sdrf_df, msstats_df)
            else:
                msstats_df = parquet_common_process(msstats_df)
            result_df = data_common_process(msstats_df, min_aa)

            # Write CSVs by Sample ID
            for sample in sample_names:
                file_name = f"{temp}/{sample}.csv"
                write_mode = "a" if os.path.exists(file_name) else "w"
                header = False if os.path.exists(file_name) else True
                result_df[result_df[SAMPLE_ID] == sample].to_csv(
                    file_name, index=False, header=header, mode=write_mode
                )
            unique_df = result_df.groupby([PEPTIDE_CANONICAL]).filter(
                lambda x: len(set(x[PROTEIN_NAME])) == 1
            )[[PEPTIDE_CANONICAL, PROTEIN_NAME]]
            unique_dict = dict(
                zip(unique_df[PEPTIDE_CANONICAL], unique_df[PROTEIN_NAME])
            )
            for i in unique_dict.keys():
                if i in unique_peptides.keys() and unique_dict[i] != unique_peptides[i]:
                    unique_peptides.pop(i)
                else:
                    unique_peptides[i] = unique_dict[i]

        proteins_list = list(unique_peptides.values())
        count_dict = {
            element: proteins_list.count(element) for element in set(proteins_list)
        }
        strong_proteins = [
            element for element in count_dict if count_dict[element] >= min_unique
        ]
        del proteins_list, count_dict
        print(f"Number of unique peptides: {len(list(unique_peptides.keys()))}")
        print(f"Number of strong proteins: {len(strong_proteins)}")

        # TODO: Filter proteins with less unique peptides than min_unique (default: 2)
        plot_samples = random.sample(sample_names, min(len(sample_names), 20))
        plot_width = 10 + len(plot_samples) * 0.5
        pdf = PdfPages(qc_report)
        original_intensities_df = pd.DataFrame()

        print("INFO: Second iteration to filter data and prepare normalization...")
        print("Logarithmic if specified..")
        norm_record = [0] * 2
        for sample in sample_names:
            msstats_df = pd.read_csv(f"{temp}/{sample}.csv", sep=",")
            msstats_df = msstats_df[msstats_df[PROTEIN_NAME].isin(strong_proteins)]
            # Remove high abundant and contaminants proteins and the outliers
            if remove_ids is not None:
                msstats_df = remove_protein_by_ids(msstats_df, remove_ids)
            if remove_decoy_contaminants:
                msstats_df = remove_contaminants_decoys(msstats_df)
                norm_record[0] += len(msstats_df)
            msstats_df.loc[msstats_df.Intensity == 0, INTENSITY] = 1
            msstats_df[NORM_INTENSITY] = (
                np.log2(msstats_df[INTENSITY]) if log2 else msstats_df[INTENSITY]
            )
            if sample in plot_samples:
                original_intensities_df = pd.concat(
                    [original_intensities_df, msstats_df]
                )
            if not skip_normalization:
                if nmethod == "msstats":
                    if label in ["TMT", "ITRAQ"]:
                        g = msstats_df.groupby(["Run", "Channel"])
                    else:
                        g = msstats_df.groupby(["Run", "Fraction"])
                    for name, group in g:
                        group_intensity = group[NORM_INTENSITY].tolist()
                        if name not in group_intensities:
                            group_intensities[name] = group_intensity
                        else:
                            group_intensities.update(
                                {
                                    name: group_intensities[NORM_INTENSITY]
                                    + group_intensity
                                }
                            )
                elif nmethod == "quantile":
                    msstats_df = (
                        msstats_df.groupby(
                            [
                                PEPTIDE_SEQUENCE,
                                PEPTIDE_CANONICAL,
                                PEPTIDE_CHARGE,
                                FRACTION,
                                RUN,
                                BIOREPLICATE,
                                PROTEIN_NAME,
                                STUDY_ID,
                                CONDITION,
                            ]
                        )[NORM_INTENSITY]
                        .agg(np.mean)
                        .reset_index()
                    )
                    rank = msstats_df[NORM_INTENSITY].rank(method="average")
                    dic = dict(zip(rank, msstats_df[NORM_INTENSITY]))
                    if len(quantile) == 0:
                        quantile = {k: (v, 1) for k, v in dic.items()}
                    else:
                        # update = min(len(quantile), len(dic))
                        intersec = set(quantile.keys()) & set(dic.keys())
                        update = set(dic.keys()) - set(quantile.keys())
                        quantile.update(
                            {
                                i: (quantile[i][0] + dic[i], quantile[i][1] + 1)
                                for i in intersec
                            }
                        )
                        if len(update) > 0:
                            quantile.update({k: (dic[k], 1) for k in update})
                    msstats_df[SAMPLE_ID] = sample
                else:
                    exit("Stream process only supports msstats and quantile methods!")
            msstats_df.to_csv(f"{temp}/{sample}.csv", index=False, sep=",")
            norm_record[1] += len(msstats_df)
        if not skip_normalization and nmethod == "quantile":
            quantile = {k: v[0] / v[1] for k, v in quantile.items()}
        print(f"Peptides after contaminants removal: {norm_record[0]}")
        print(f"Number of peptides after normalization: {norm_record[1]}")
        # Save original intensities QC plots
        original_intensities_df = original_intensities_df.reset_index(drop=True)
        density = plot_distributions(
            original_intensities_df,
            INTENSITY,
            SAMPLE_ID,
            log2=not log2,
            width=plot_width,
            title="Original peptidoform intensity distribution (no normalization)",
        )
        pdf.savefig(density)
        box = plot_box_plot(
            original_intensities_df,
            INTENSITY,
            SAMPLE_ID,
            log2=not log2,
            width=plot_width,
            title="Original peptidoform intensity distribution (no normalization)",
            violin=violin,
        )
        plt.show()
        pdf.savefig(box)
        del original_intensities_df

        # TODO: Peptide intensity normalization
        peptides_count = pd.DataFrame(
            columns=[PROTEIN_NAME, PEPTIDE_CANONICAL, "count"]
        )
        norm_intensities_df = pd.DataFrame()
        if not skip_normalization and nmethod == "msstats":
            # For ISO normalization
            if label in ["TMT", "ITRAQ"]:
                median_baseline = np.median(
                    list(set(sum(group_intensities.values(), [])))
                )
                group_intensities = {
                    key: np.median(list(values))
                    for key, values in group_intensities.items()
                }
            else:
                fractions = [i[1] for i in group_intensities.keys()]
                fraction_median = {}
                for fraction in fractions:
                    fraction_keys = [
                        i for i in group_intensities.keys() if i[1] == fraction
                    ]
                    fraction_intensities = []
                    for key in fraction_keys:
                        fraction_intensities.extend(group_intensities[key])
                    fraction_median[fraction] = np.median(fraction_intensities)
                group_intensities = {
                    key: np.median(values) for key, values in group_intensities.items()
                }
        print("INFO: Third iteration to normalize and counting peptides frequency...")
        size_record = [0] * 3

        def normalization(
            dataset_df, label, sample, skip_normalization, nmethod, record
        ):
            if not skip_normalization:
                field = NORM_INTENSITY
                if nmethod == "msstats":
                    # For ISO normalization
                    if label in ["TMT", "ITRAQ"]:
                        dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(
                            lambda x: x[field]
                            - group_intensities[(x["Run"], x["Channel"])]
                            + median_baseline,
                            axis=1,
                        )
                    else:
                        dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(
                            lambda x: x[field]
                            - group_intensities[(x["Run"], x["Fraction"])]
                            + np.median(
                                [
                                    group_intensities[i]
                                    for i in group_intensities.keys()
                                    if i[1] == x["Fraction"]
                                ]
                            ),
                            axis=1,
                        )
                elif nmethod == "quantile":
                    rank = dataset_df[NORM_INTENSITY].rank(method="average")
                    ref_dict = dict(zip(rank, dataset_df[NORM_INTENSITY]))
                    ref_dict = {v: quantile[k] for k, v in ref_dict.items()}
                    dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(
                        lambda x: ref_dict.get(x[NORM_INTENSITY], np.nan),
                        axis=1,
                    )
            dataset_df = dataset_df.drop_duplicates()
            dataset_df = dataset_df[dataset_df[NORM_INTENSITY].notna()]
            dataset_df = get_peptidoform_normalize_intensities(dataset_df)
            record[0] += len(dataset_df.index)
            dataset_df = sum_peptidoform_intensities(dataset_df)
            record[1] += len(dataset_df.index)
            dataset_df = average_peptide_intensities(dataset_df)
            record[2] += len(dataset_df.index)

            return dataset_df, record

        for sample in sample_names:
            dataset_df = pd.read_csv(f"{temp}/{sample}.csv", sep=",")
            if len(dataset_df) != 0:
                norm_df, size_record = normalization(
                    dataset_df, label, sample, skip_normalization, nmethod, size_record
                )
            else:
                continue
            sample_peptides = norm_df[PEPTIDE_CANONICAL].unique().tolist()
            if remove_low_frequency_peptides and sample_number > 1:
                sample_peptides = norm_df[
                    [PROTEIN_NAME, PEPTIDE_CANONICAL]
                ].drop_duplicates()
                sample_peptides["count"] = 1
                peptides_count = (
                    pd.concat([peptides_count, sample_peptides])
                    .groupby([PROTEIN_NAME, PEPTIDE_CANONICAL])
                    .agg(sum)
                    .reset_index()
                )
            norm_df.to_csv(f"{temp}/{sample}.csv", sep=",", index=False)
            if sample in plot_samples:
                norm_intensities_df = pd.concat([norm_intensities_df, norm_df])
        del group_intensities, quantile
        print(f"Number of peptides after peptidofrom selection: {size_record[0]}")
        print(f"Number of peptides after selection: {size_record[1]}")
        print(f"Number of peptides after average: {size_record[2]}")
        # Save normalized intensities QC plots
        norm_intensities_df = norm_intensities_df.reset_index(drop=True)
        log_after_norm = nmethod == "msstats" or (
            (nmethod == "quantile" or nmethod == "robust") and not log2
        )
        density = plot_distributions(
            norm_intensities_df,
            NORM_INTENSITY,
            SAMPLE_ID,
            log2=log_after_norm,
            width=plot_width,
            title="Peptidoform intensity distribution after normalization, method: "
            + nmethod,
        )
        plt.show()
        pdf.savefig(density)
        box = plot_box_plot(
            norm_intensities_df,
            NORM_INTENSITY,
            SAMPLE_ID,
            log2=log_after_norm,
            width=plot_width,
            title="Peptidoform intensity distribution after normalization, method: "
            + nmethod,
            violin=violin,
        )
        plt.show()
        pdf.savefig(box)
        del norm_intensities_df, strong_proteins

        print("INFO: Writing normalized intensities into CSV...")
        if remove_low_frequency_peptides and sample_number > 1:
            peptides_count = peptides_count.loc[
                (peptides_count["count"] > 0.20 * sample_number)
                & (peptides_count["count"] != sample_number - 1)
            ]

        final_norm_intensities_df = pd.DataFrame()
        size_record = 0
        for sample in sample_names:
            dataset_df = pd.read_csv(f"{temp}/{sample}.csv", sep=",")
            if remove_low_frequency_peptides and sample_number > 1:
                # Filter low-frequency peptides, which indicate whether the peptide occurs less than 20% in all samples or
                # only in one sample
                dataset_df = dataset_df.merge(
                    peptides_count[[PEPTIDE_CANONICAL, PROTEIN_NAME]], how="inner"
                )
            size_record += len(dataset_df.index)
            dataset_df = dataset_df[
                [PEPTIDE_CANONICAL, PROTEIN_NAME, SAMPLE_ID, NORM_INTENSITY, CONDITION]
            ]
            write_mode = "a" if os.path.exists(output) else "w"
            header = False if os.path.exists(output) else True
            dataset_df.to_csv(output, index=False, header=header, mode=write_mode)
            dataset_df.to_csv(f"{temp}/{sample}.csv", sep=",", index=False)
            if sample in plot_samples:
                final_norm_intensities_df = pd.concat(
                    [final_norm_intensities_df, dataset_df]
                )
        print(f"Peptides after remove low frequency peptides: {size_record}")
        if remove_low_frequency_peptides:
            del peptides_count

        # TODO: No peptides intensity normalization applied in stream processing.
        # Save final normalized intensities QC plots
        log_after_norm = nmethod == "msstats" or (
            (nmethod == "quantile" or nmethod == "robust") and not log2
        )
        final_norm_intensities_df = final_norm_intensities_df.reset_index(drop=True)
        density = plot_distributions(
            final_norm_intensities_df,
            NORM_INTENSITY,
            SAMPLE_ID,
            log2=log_after_norm,
            width=plot_width,
            title="Normalization at peptide level method: " + nmethod,
        )
        plt.show()
        pdf.savefig(density)
        box = plot_box_plot(
            final_norm_intensities_df,
            NORM_INTENSITY,
            SAMPLE_ID,
            log2=log_after_norm,
            width=plot_width,
            title="Normalization at peptide level method: " + nmethod,
            violin=violin,
        )
        plt.show()
        pdf.savefig(box)
        pdf.close()


if __name__ == "__main__":
    peptide_normalization()
