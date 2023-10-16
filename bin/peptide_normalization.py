#!/usr/bin/env python
import gc

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qnorm
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame

from ibaq.ibaqpy_commons import (
    BIOREPLICATE,
    CHANNEL,
    CONDITION,
    FEATURE_COLUMNS,
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
    ITRAQ4plex,
    ITRAQ8plex,
    TMT6plex,
    TMT10plex,
    TMT11plex,
    TMT16plex,
    get_canonical_peptide,
    get_reference_name,
    get_study_accession,
    parquet_map,
    parse_uniprot_accession,
    plot_box_plot,
    plot_distributions,
    remove_contaminants_decoys,
    remove_extension_file,
    remove_protein_by_ids,
    sum_peptidoform_intensities,
)


def print_dataset_size(dataset: DataFrame, message: str, verbose: bool) -> None:
    if verbose:
        print(message + str(len(dataset.index)))


def print_help_msg(command) -> None:
    """
    Print help information
    :param command: command to print helps
    :return: print help
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


def remove_outliers_iqr(dataset: DataFrame):
    """
    This method removes outliers from the dataframe inplace, the variable used for the outlier removal is Intensity
    :param dataset: Peptide dataframe
    :return: None
    """
    q1 = dataset[INTENSITY].quantile(0.25)
    q3 = dataset[INTENSITY].quantile(0.75)
    iqr = q3 - q1

    dataset.query("(@q1 - 1.5 * @iqr) <= Intensity <= (@q3 + 1.5 * @iqr)", inplace=True)


def remove_missing_values(normalize_df: DataFrame, ratio: float = 0.3) -> DataFrame:
    """
    Remove missing values if the peptide do not have values in the most of the samples
    :param normalize_df: data frame with the data
    :param ratio: ratio of samples without intensity values.
    :return:
    """
    n_samples = len(normalize_df.columns)
    normalize_df = normalize_df.dropna(thresh=round(n_samples * ratio))
    return normalize_df


def intensity_normalization(
    dataset: DataFrame,
    field: str,
    class_field: str = "all",
    scaling_method: str = "msstats",
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

    elif scaling_method == "qnorm":
        # pivot to have one col per sample
        print("Transforming to wide format dataset size {}".format(len(dataset.index)))
        normalize_df = pd.pivot_table(
            dataset,
            index=[
                PEPTIDE_SEQUENCE,
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
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(
            id_vars=[
                PEPTIDE_SEQUENCE,
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
        print(dataset.head())
        return normalize_df

    return dataset


def get_peptidoform_normalize_intensities(
    dataset: DataFrame, higher_intensity: bool = True
) -> DataFrame:
    """
    Select the best peptidoform for the same sample and the same replicates. A peptidoform is the combination of
    a (PeptideSequence + Modifications) + Charge state.
    :param dataset: dataset including all properties
    :param higher_intensity: select based on normalize intensity, if false based on best scored peptide
    :return:
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    if higher_intensity:
        dataset = dataset.loc[
            dataset.groupby(
                [PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE],
                observed=True,
            )[NORM_INTENSITY].idxmax()
        ].reset_index(drop=True)
    else:
        dataset = dataset.loc[
            dataset.groupby(
                [PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE],
                observed=True,
            )[SEARCH_ENGINE].idxmax()
        ].reset_index(drop=True)
    return dataset


def average_peptide_intensities(dataset: DataFrame) -> DataFrame:
    """
    Median the intensities of all the peptidoforms for a specific peptide sample combination.
    :param dataset: Dataframe containing all the peptidoforms
    :return: New dataframe
    """
    dataset_df = dataset.groupby(
        [PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION], observed=True
    )[NORM_INTENSITY].median()
    dataset_df = dataset_df.reset_index()
    dataset_df = pd.merge(
        dataset_df,
        dataset[[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION]],
        how="left",
        on=[PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION],
    )
    return dataset_df


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
    if scaling_method == "qnorm":
        # pivot to have one col per sample
        normalize_df = pd.pivot_table(
            dataset_df,
            index=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION],
            columns=class_field,
            values=field,
            aggfunc={field: np.mean},
            observed=True,
        )
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(
            id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION]
        )
        normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
        normalize_df = normalize_df[normalize_df[NORM_INTENSITY].notna()]
        return normalize_df

    return dataset_df


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
    help="Normalization method used to normalize intensities for all samples (options: qnorm)",
    default="qnorm",
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

    compression_method = "gzip" if compress else None

    if parquet is None:
        # Read the msstats file
        feature_df = pd.read_csv(
            msstats,
            sep=",",
            compression=compression_method,
            dtype={CONDITION: "category", ISOTOPE_LABEL_TYPE: "category"},
        )

        feature_df.rename(
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

        feature_df[PROTEIN_NAME] = feature_df[PROTEIN_NAME].apply(
            parse_uniprot_accession
        )

        # Read the sdrf file
        sdrf_df = pd.read_csv(sdrf, sep="\t", compression=compression_method)
        sdrf_df[REFERENCE] = sdrf_df["comment[data file]"].apply(remove_extension_file)
        print(sdrf_df)

        if FRACTION not in feature_df.columns:
            feature_df[FRACTION] = 1
            feature_df = feature_df[
                [
                    PROTEIN_NAME,
                    PEPTIDE_SEQUENCE,
                    PEPTIDE_CHARGE,
                    INTENSITY,
                    REFERENCE,
                    CONDITION,
                    RUN,
                    BIOREPLICATE,
                    FRACTION,
                    FRAGMENT_ION,
                    ISOTOPE_LABEL_TYPE,
                ]
            ]

        # Merged the SDRF with the Resulted file
        labels = set(sdrf_df["comment[label]"])
        if CHANNEL not in feature_df.columns:
            feature_df[REFERENCE] = feature_df[REFERENCE].apply(remove_extension_file)
            dataset_df = pd.merge(
                feature_df,
                sdrf_df[["source name", REFERENCE]],
                how="left",
                on=[REFERENCE],
            )
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
            choice = (
                pd.DataFrame.from_dict(choice, orient="index", columns=[CHANNEL])
                .reset_index()
                .rename(columns={"index": "comment[label]"})
            )
            sdrf_df = sdrf_df.merge(choice, on="comment[label]", how="left")
            feature_df[REFERENCE] = feature_df[REFERENCE].apply(get_reference_name)
            dataset_df = pd.merge(
                feature_df,
                sdrf_df[["source name", REFERENCE, CHANNEL]],
                how="left",
                on=[REFERENCE, CHANNEL],
            )
            # result_df.drop(CHANNEL, axis=1, inplace=True)
            dataset_df = dataset_df[dataset_df["Condition"] != "Empty"]
            dataset_df.rename(columns={"Charge": PEPTIDE_CHARGE}, inplace=True)
        elif "ITRAQ" in ",".join(labels) or "itraq" in ",".join(labels):
            if len(labels) > 4:
                choice = ITRAQ8plex
            else:
                choice = ITRAQ4plex
            choice = (
                pd.DataFrame.from_dict(choice, orient="index", columns=[CHANNEL])
                .reset_index()
                .rename(columns={"index": "comment[label]"})
            )
            sdrf_df = sdrf_df.merge(choice, on="comment[label]", how="left")
            feature_df[REFERENCE] = feature_df[REFERENCE].apply(get_reference_name)
            dataset_df = pd.merge(
                feature_df,
                sdrf_df[["source name", REFERENCE, CHANNEL]],
                how="left",
                on=[REFERENCE, CHANNEL],
            )
            dataset_df = dataset_df[dataset_df["Condition"] != "Empty"]
            dataset_df.rename(columns={"Charge": PEPTIDE_CHARGE}, inplace=True)
        else:
            print("Warning: Only support label free, TMT and ITRAQ experiment!")
            exit(1)

        # Remove the intermediate variables and free the memory
        del feature_df, sdrf_df
        gc.collect()
    else:
        dataset_df = pd.read_parquet(parquet)[FEATURE_COLUMNS]
        dataset_df = dataset_df.rename(columns=parquet_map)
        dataset_df[PROTEIN_NAME] = dataset_df.apply(
            lambda x: ",".join(x[PROTEIN_NAME]), axis=1
        )
        label_type = dataset_df[CHANNEL].unique().tolist()
        if len(label_type) == 1:
            dataset_df.drop(CHANNEL, inplace=True, axis=1)
        dataset_df = dataset_df[dataset_df["Condition"] != "Empty"]

    # Remove 0 intensity signals from the msstats file
    dataset_df = dataset_df[dataset_df[INTENSITY] > 0]
    dataset_df[PEPTIDE_CANONICAL] = dataset_df.apply(
        lambda x: get_canonical_peptide(x[PEPTIDE_SEQUENCE]), axis=1
    )
    # Only peptides with more than min_aa (default: 7) amino acids are retained
    dataset_df = dataset_df[
        dataset_df.apply(lambda x: len(x[PEPTIDE_CANONICAL]) >= min_aa, axis=1)
    ]
    # Only proteins with unique peptides number greater than min_unique (default: 2) are retained
    unique_peptides = set(
        dataset_df.groupby(PEPTIDE_CANONICAL)
        .filter(lambda x: len(set(x[PROTEIN_NAME])) == 1)[PEPTIDE_CANONICAL]
        .tolist()
    )
    strong_proteins = set(
        dataset_df[dataset_df[PEPTIDE_CANONICAL].isin(unique_peptides)]
        .groupby(PROTEIN_NAME)
        .filter(lambda x: len(set(x[PEPTIDE_CANONICAL])) >= min_unique)[PROTEIN_NAME]
        .tolist()
    )
    dataset_df = dataset_df[dataset_df[PROTEIN_NAME].isin(strong_proteins)]

    if msstats:
        dataset_df.rename(columns={"source name": SAMPLE_ID}, inplace=True)
    dataset_df[STUDY_ID] = dataset_df[SAMPLE_ID].apply(get_study_accession)
    dataset_df = dataset_df.filter(
        items=[
            PEPTIDE_SEQUENCE,
            PEPTIDE_CHARGE,
            FRACTION,
            RUN,
            BIOREPLICATE,
            PROTEIN_NAME,
            STUDY_ID,
            CONDITION,
            SAMPLE_ID,
            INTENSITY,
        ]
    )
    dataset_df[CONDITION] = pd.Categorical(dataset_df[CONDITION])
    dataset_df[STUDY_ID] = pd.Categorical(dataset_df[STUDY_ID])
    dataset_df[SAMPLE_ID] = pd.Categorical(dataset_df[SAMPLE_ID])

    pd.set_option("display.max_columns", None)
    print("Loading data..")
    print_dataset_size(dataset_df, "Number of peptides: ", verbose)

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

    print("Select the best peptidoform across fractions...")
    print(
        "Number of peptides before peptidofrom selection: " + str(len(dataset_df.index))
    )
    dataset_df = get_peptidoform_normalize_intensities(dataset_df)
    print(
        "Number of peptides after peptidofrom selection: " + str(len(dataset_df.index))
    )

    # Add the peptide sequence canonical without the modifications
    if PEPTIDE_CANONICAL not in dataset_df.columns:
        print("Add Canonical peptides to the dataframe...")
        dataset_df[PEPTIDE_CANONICAL] = dataset_df[PEPTIDE_SEQUENCE].apply(
            lambda x: get_canonical_peptide(x)
        )

    print("Sum all peptidoforms per Sample...")
    print("Number of peptides before sum selection: " + str(len(dataset_df.index)))
    dataset_df = sum_peptidoform_intensities(dataset_df)
    print("Number of peptides after sum: " + str(len(dataset_df.index)))

    print("Average all peptidoforms per Peptide/Sample...")
    print("Number of peptides before average: " + str(len(dataset_df.index)))
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

    if remove_low_frequency_peptides:
        print(
            "Peptides before removing low frequency peptides: "
            + str(len(dataset_df.index))
        )
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


if __name__ == "__main__":
    peptide_normalization()
