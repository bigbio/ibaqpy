import os
import re

import click
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pandas import DataFrame

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


def remove_protein_by_ids(
    dataset: DataFrame, protein_file: str, protein_field=PROTEIN_NAME
) -> DataFrame:
    """
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param protein_file: contaminants file
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants_reader = open(protein_file, "r")
    contaminants = contaminants_reader.read().split("\n")
    contaminants = [cont for cont in contaminants if cont.strip()]
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex)]


def remove_contaminants_decoys(
    dataset: DataFrame, protein_field=PROTEIN_NAME
) -> DataFrame:
    """
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param contaminants_file: contaminants file
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants = []
    contaminants.append("CONTAMINANT")
    contaminants.append("DECOY")
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex)]


def get_canonical_peptide(peptide_sequence: str) -> str:
    """
    This function returns a peptide sequence without the modification information
    :param peptide_sequence: peptide sequence with mods
    :return: peptide sequence
    """
    clean_peptide = re.sub("[\(\[].*?[\)\]]", "", peptide_sequence)
    clean_peptide = clean_peptide.replace(".", "").replace("-", "")
    return clean_peptide


def plot_distributions(
    dataset: DataFrame,
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
    dataset: DataFrame,
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


def sum_peptidoform_intensities(dataset: DataFrame) -> DataFrame:
    """
    Sum the peptidoform intensities for all peptidofrom across replicates of the same sample.
    :param dataset: Dataframe to be analyzed
    :return: dataframe with the intensities
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    normalize_df = dataset.groupby(
        [PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION], observed=True
    )[NORM_INTENSITY].sum()
    normalize_df = normalize_df.reset_index()
    normalize_df = pd.merge(
        normalize_df,
        dataset[[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION]],
        how="left",
        on=[PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION],
    )
    normalize_df.drop_duplicates(inplace=True)
    return normalize_df


def parse_uniprot_accession(uniprot_id: str) -> str:
    """
    Parse the uniprot accession from the uniprot id in the form of
    tr|CONTAMINANT_Q3SX28|CONTAMINANT_TPM2_BOVIN and convert to CONTAMINANT_Q3SX28
    :param uniprot_id: uniprot id
    :return: uniprot accession
    """
    uniprot_list = uniprot_id.split(";")
    result_uniprot_list = []
    for accession in uniprot_list:
        if accession.count("|") == 2:
            accession = accession.split("|")[1]
        result_uniprot_list.append(accession)
    return ";".join(result_uniprot_list)


def get_study_accession(sample_id: str) -> str:
    """
    Get the project accession from the Sample accession. The function expected a sample accession in the following
    format PROJECT-SAMPLEID
    :param sample_id: Sample Accession
    :return: study accession
    """
    return sample_id.split("-")[0]


def get_spectrum_prefix(reference_spectrum: str) -> str:
    """
    Get the reference name from Reference column. The function expected a reference name in the following format eg.
    20150820_Haura-Pilot-TMT1-bRPLC03-2.mzML_controllerType=0 controllerNumber=1 scan=16340. This function can also
    remove suffix of spectrum files.
    :param reference_spectrum: 
    :return: reference name
    """
    return re.split(r"\.mzML|\.MZML|\.raw|\.RAW|\.d|\.wiff", reference_spectrum)[0]


# Common functions when normalizing peptide dataframe
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
    dataset_df.drop_duplicates(inplace=True)
    return dataset_df


# Functions needed by Combiner
def load_sdrf(sdrf_path: str) -> DataFrame:
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


def load_feature(feature_path: str) -> DataFrame:
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
