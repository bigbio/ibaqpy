import re
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from pandas import DataFrame

PROTEIN_NAME = 'ProteinName'
PEPTIDE_SEQUENCE = 'PeptideSequence'
PEPTIDE_CANONICAL = "PeptideCanonical"
PEPTIDE_CHARGE = 'PrecursorCharge'
FRAGMENT_ION = 'FragmentIon'
PRODUCT_CHARGE = 'ProductCharge'
ISOTOPE_LABEL_TYPE = 'IsotopeLabelType'
CHANNEL = 'Channel'
MIXTRUE = 'Mixture'
TECHREPMIXTURE = 'TechRepMixture'
CONDITION = 'Condition'
BIOREPLICATE = 'BioReplicate'
RUN = 'Run'
FRACTION = 'Fraction'
INTENSITY = 'Intensity'
NORM_INTENSITY = 'NormIntensity'
RT = 'Rt'
REFERENCE = 'Reference'
SAMPLE_ID = 'SampleID'
STUDY_ID = 'StudyID'
SEARCH_ENGINE = 'searchScore'
SCAN = 'Scan'
MBR = 'MatchBetweenRuns'
IBAQ = 'Ibaq'
IBAQ_NORMALIZED = 'IbaqNorm'
IBAQ_LOG = 'IbaqLog'
IBAQ_PPB = 'IbaqPpb'

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

TMT6plex = {"TMT126": 1, "TMT127": 2, "TMT128": 3, "TMT129": 4, "TMT130": 5, "TMT131": 6}

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


def remove_contaminants_decoys(dataset: DataFrame, contaminants_file: str, protein_field=PROTEIN_NAME) -> DataFrame:
    """
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param contaminants_file: contaminants file
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants_reader = open(contaminants_file, 'r')
    contaminants = contaminants_reader.read().split("\n")
    contaminants = [cont for cont in contaminants if cont.strip()]

    contaminants.append('CONTAMINANT')
    contaminants.append('DECOY')
    # cregex = ".*(" + '|'.join(contaminants) + ").*"
    cregex = '|'.join(contaminants)
    # for contaminant in contaminants:
    # dataset.drop(index=dataset[dataset[protein_field].str.contains(contaminant)].index, inplace=True)

    return dataset[~dataset[protein_field].str.contains(cregex)]


def get_canonical_peptide(peptide_sequence: str) -> str:
    """
    This function returns a peptide sequence without the modification information
    :param peptide_sequence: peptide sequence with mods
    :return: peptide sequence
    """
    clean_peptide = re.sub("[\(\[].*?[\)\]]", "", peptide_sequence)
    clean_peptide = clean_peptide.replace(".", "")
    return clean_peptide


def plot_distributions(dataset: DataFrame, field: str, class_field: str, log2: bool = True) -> None:
    """
    Print the quantile plot for the dataset
    :param dataset: DataFrame
    :param field: Field that would be use in the dataframe to plot the quantile
    :param class_field: Field to group the quantile into classes
    :param log2: Log the intensity values
    :return:
    """
    pd.set_option('mode.chained_assignment', None)
    normalize = dataset[[field, class_field]]
    if log2:
        normalize[field] = np.log2(normalize[field])
    normalize.dropna(subset=[field], inplace=True)
    data_wide = normalize.pivot(columns=class_field,
                                values=field)
    # plotting multiple density plot
    data_wide.plot.kde(figsize=(8, 6), linewidth=2, legend=False)
    pd.set_option('mode.chained_assignment', 'warn')


def plot_box_plot(dataset: DataFrame, field: str, class_field: str, log2: bool = False, weigth: int = 10,
                  rotation: int = 45, title: str = "", violin: bool = False) -> None:
    """
    Plot a box plot of two values field and classes field
    :param violin: Also add violin on top of box plot
    :param dataset: Dataframe with peptide intensities
    :param field: Intensity field
    :param class_field: class to group the peptides
    :param log2: transform peptide intensities to log scale
    :param weigth: size of the plot
    :param rotation: rotation of the x-axis
    :param title: Title of the box plot
    :return:
    """
    pd.set_option('mode.chained_assignment', None)
    normalized = dataset[[field, class_field]]
    np.seterr(divide='ignore')
    plt.figure(figsize=(weigth, 10))
    if log2:
        normalized[field] = np.log2(normalized[field])

    if violin:
        chart = sns.violinplot(x=class_field, y=field, data=normalized, boxprops=dict(alpha=.3), palette="muted")
    else:
        chart = sns.boxplot(x=class_field, y=field, data=normalized, boxprops=dict(alpha=.3), palette="muted")

    chart.set(title=title)
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation)
    plt.show()
    pd.set_option('mode.chained_assignment', 'warn')
