#!/usr/bin/env python
import os
import click
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

PARQUET_COLUMNS = [
    "protein_accessions",
    "peptidoform",
    "sequence",
    "charge",
    "fragment_ion",
    "isotope_label_type",
    "channel",
    "condition",
    "biological_replicate",
    "run",
    "fraction",
    "intensity",
    "reference_file_name",
    "sample_accession",
]

PROTEIN_NAME = "ProteinName"
PEPTIDE_SEQUENCE = "PeptideSequence"
PEPTIDE_CANONICAL = "PeptideCanonical"
PEPTIDE_CHARGE = "PrecursorCharge"
FRAGMENT_ION = "FragmentIon"
PRODUCT_CHARGE = "ProductCharge"
ISOTOPE_LABEL_TYPE = "IsotopeLabelType"
CHANNEL = "Channel"
MIXTRUE = "Mixture"
TECHREPMIXTURE = "TechRepMixture"
CONDITION = "Condition"
BIOREPLICATE = "BioReplicate"
TECHREPLICATE = "TechReplicate"
RUN = "Run"
FRACTION = "Fraction"
INTENSITY = "Intensity"
NORM_INTENSITY = "NormIntensity"
RT = "Rt"
REFERENCE = "Reference"
SAMPLE_ID = "SampleID"
STUDY_ID = "StudyID"
SEARCH_ENGINE = "searchScore"
SCAN = "Scan"
MBR = "MatchBetweenRuns"
IBAQ = "Ibaq"
IBAQ_NORMALIZED = "IbaqNorm"
IBAQ_LOG = "IbaqLog"
IBAQ_PPB = "IbaqPpb"

parquet_map = {
    "protein_accessions": PROTEIN_NAME,
    "peptidoform": PEPTIDE_SEQUENCE,
    "sequence": PEPTIDE_CANONICAL,
    "charge": PEPTIDE_CHARGE,
    "fragment_ion": FRAGMENT_ION,
    "isotope_label_type": ISOTOPE_LABEL_TYPE,
    "channel": CHANNEL,
    "condition": CONDITION,
    "biological_replicate": BIOREPLICATE,
    "run": RUN,
    "fraction": FRACTION,
    "intensity": INTENSITY,
    "reference_file_name": REFERENCE,
    "sample_accession": SAMPLE_ID,
}

TMT16plex = {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131N": 10,
    "TMT131C": 11,
    "TMT132N": 12,
    "TMT132C": 13,
    "TMT133N": 14,
    "TMT133C": 15,
    "TMT134N": 16,
}

TMT11plex = {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131N": 10,
    "TMT131C": 11,
}

TMT10plex = {
    "TMT126": 1,
    "TMT127N": 2,
    "TMT127C": 3,
    "TMT128N": 4,
    "TMT128C": 5,
    "TMT129N": 6,
    "TMT129C": 7,
    "TMT130N": 8,
    "TMT130C": 9,
    "TMT131": 10,
}

TMT6plex = {
    "TMT126": 1,
    "TMT127": 2,
    "TMT128": 3,
    "TMT129": 4,
    "TMT130": 5,
    "TMT131": 6,
}

ITRAQ4plex = {"ITRAQ114": 1, "ITRAQ115": 2, "ITRAQ116": 3, "ITRAQ117": 4}

ITRAQ8plex = {
    "ITRAQ113": 1,
    "ITRAQ114": 2,
    "ITRAQ115": 3,
    "ITRAQ116": 4,
    "ITRAQ117": 5,
    "ITRAQ118": 6,
    "ITRAQ119": 7,
    "ITRAQ121": 8,
}


def print_help_msg(command: click.Command):
    """
    Print the help of the command
    :param command: click command object
    :return: None
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


def get_accession(identifier: str) -> str:
    """
    Get protein accession from the identifier  (e.g. sp|P12345|PROT_NAME)
    :param identifier: Protein identifier
    :return: Protein accession
    """
    identifier_lst = identifier.split("|")
    if len(identifier_lst) == 1:
        return identifier_lst[0]
    else:
        return identifier_lst[1]


def plot_distributions(
    dataset: pd.DataFrame,
    field: str,
    class_field: str,
    title: str = "",
    log2: bool = True,
    width: float = 10,
) -> matplotlib.pyplot:
    """
    Print the quantile plot for the dataset
    :param dataset: DataFrame
    :param field: Field that would be use in the dataframe to plot the quantile
    :param class_field: Field to group the quantile into classes
    :param title: Title of the box plot
    :param log2: Log the intensity values
    :param width: size of the plot
    :return:
    """
    pd.set_option("mode.chained_assignment", None)
    normalize = dataset[[field, class_field]].reset_index(drop=True)
    if log2:
        normalize[field] = np.log2(normalize[field])
    normalize.dropna(subset=[field], inplace=True)
    data_wide = normalize.pivot(columns=class_field, values=field)
    # plotting multiple density plot
    data_wide.plot.kde(figsize=(width, 8), linewidth=2, legend=False)
    plt.title(title)
    pd.set_option("mode.chained_assignment", "warn")

    return plt.gcf()


def plot_box_plot(
    dataset: pd.DataFrame,
    field: str,
    class_field: str,
    log2: bool = False,
    width: float = 10,
    rotation: int = 30,
    title: str = "",
    violin: bool = False,
) -> matplotlib.pyplot:
    """
    Plot a box plot of two values field and classes field
    :param violin: Also add violin on top of box plot
    :param dataset: Dataframe with peptide intensities
    :param field: Intensity field
    :param class_field: class to group the peptides
    :param log2: transform peptide intensities to log scale
    :param width: size of the plot
    :param rotation: rotation of the x-axis
    :param title: Title of the box plot
    :return:
    """
    pd.set_option("mode.chained_assignment", None)
    normalized = dataset[[field, class_field]]
    np.seterr(divide="ignore")
    plt.figure(figsize=(width, 14))
    if log2:
        normalized[field] = np.log2(normalized[field])

    if violin:
        chart = sns.violinplot(
            x=class_field,
            y=field,
            data=normalized,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )
    else:
        chart = sns.boxplot(
            x=class_field,
            y=field,
            data=normalized,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )

    chart.set(title=title)
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation, ha="right")
    pd.set_option("mode.chained_assignment", "warn")

    return plt.gcf()


# Functions needed by Combiner
def load_sdrf(sdrf_path: str) -> pd.DataFrame:
    """
    Load sdrf TSV as a dataframe.
    :param sdrf_path: Path to SDRF TSV.
    :return:
    """
    if not os.path.exists(sdrf_path):
        raise FileNotFoundError(f"{sdrf_path} does not exist!")
    sdrf_df = pd.read_csv(sdrf_path, sep="\t")
    sdrf_df.columns = [col.lower() for col in sdrf_df.columns]
    return sdrf_df


def load_feature(feature_path: str) -> pd.DataFrame:
    """
    Load feature file as a dataframe.
    :param feature_path: Path to feature file.
    :return:
    """
    suffix = os.path.splitext(feature_path)[1][1:]
    if suffix == "parquet":
        return pd.read_parquet(feature_path)
    elif suffix == "csv":
        return pd.read_csv(feature_path)
    else:
        raise SystemExit(
            f"{suffix} is not allowed as input, please provide msstats_in or feature parquet."
        )


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
            normalize_df = normalize_df.append(group_normalize_df, ignore_index=True)
        # else:
        #     print("nothing")
        #     # Impute the missing values
        #     imputer = MissForest(max_iter=5)
        #     imputed_data = imputer.fit_transform(group_normalize_df)
        #     group_normalize_df = pd.DataFrame(
        #         imputed_data,
        #         columns=group_normalize_df.columns,
        #         index=group_normalize_df.index,
        #     )
        #     # Melt the dataframe
        #     group_normalize_df = group_normalize_df.reset_index()
        #     group_normalize_df = group_normalize_df.melt(
        #         id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION]
        #     )
        #     group_normalize_df.rename(columns={"value": NORM_INTENSITY}, inplace=True)
        #     normalize_df = normalize_df.append(group_normalize_df, ignore_index=True)

    return normalize_df
