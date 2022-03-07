import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from seaborn_qqplot import pplot
from pandas import DataFrame
from sklearn.preprocessing import MinMaxScaler, StandardScaler, MaxAbsScaler, QuantileTransformer
from sklearn.preprocessing import RobustScaler
from scipy.stats import gamma
from sklearn.impute import KNNImputer

from ibaqpy_commons import remove_contaminants_decoys, INTENSITY, REFERENCE, SAMPLE_ID, NORM_INTENSITY, \
    PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, CONDITION
import click

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


def remove_outliers_iqr(dataset: DataFrame) -> DataFrame:
    """
    This method removes outliers from the dataframe, the variable use for the outliers removal is Intesity
    :param dataset: Peptide dataframe
    :return:
    """
    Q1 = dataset[INTENSITY].quantile(0.25)
    Q3 = dataset[INTENSITY].quantile(0.75)
    IQR = Q3 - Q1

    dataset = dataset.query('(@Q1 - 1.5 * @IQR) <= Intensity <= (@Q3 + 1.5 * @IQR)')
    return dataset


def plot_distributions(dataset: DataFrame, field: str, class_field: str, log2: bool = True, weigth : int = 10) -> None:
    """
    Print the quantile plot for the dataset
    :param dataset: DataFrame
    :param field: Field that would be use in the dataframe to plot the quantile
    :param class_field: Field to group the quantile into classes
    :param weigth: size of the plot
    :return:
    """
    normalize = dataset[[field, class_field]]
    if log2:
        normalize[field] = np.log2(normalize[field])
    data_wide = normalize.pivot(columns=class_field,
                         values=field)
    # plotting multiple density plot
    data_wide.plot.kde(figsize=(8, 6),linewidth=4, legend=False)

def plot_box_plot(dataset: DataFrame, field: str, class_field: str, log2: bool = False, weigth: int = 10, rotation: int = 45) -> None:
    """
    Plot a box plot of two values field and classes field
    :param dataset: Dataframe with peptide intensities
    :param field: Intensity field
    :param class_field: class to group the peptides
    :param log2: transform peptide intensities to log scale
    :param weigth: size of the plot
    :param rotation: rotation of the x-axis
    :return:
    """
    normalized = dataset[[field, class_field]]
    np.seterr(divide='ignore')
    plt.figure(figsize=(weigth, 10))
    if log2:
        normalized[field] = np.log2(normalized[field])

    chart = sns.boxplot(x=class_field, y=field, data=normalized, boxprops=dict(alpha=.3), palette="muted")
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation)
    plt.show()


def intensity_normalization(dataset: DataFrame, field: str, class_field: str = "all", method: str = "quantile", imputation: bool = False,
                            remove_peptides: bool = False) -> DataFrame:

    scaler = QuantileTransformer(output_distribution="normal")
    if method is "robusts":
        scaler = RobustScaler()
    if method is "minmax":
        scaler = MinMaxScaler()
    if method is "standard":
        scaler = StandardScaler()
    if method is 'maxabs':
        scaler = MaxAbsScaler()

    if class_field is "all":
        normalize_df = dataset[[field]]
        dataset[NORM_INTENSITY] = scaler.fit_transform(normalize_df)
    else:
        normalize_df = dataset[[PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, CONDITION, INTENSITY, class_field]]
        normalize_df = pd.pivot_table(normalize_df, values=[INTENSITY], index=[PEPTIDE_CHARGE, PEPTIDE_SEQUENCE, CONDITION], columns=[class_field],
                       aggfunc={INTENSITY: np.mean})

        # Remove peptides that are not found in 50% of the samples
        if remove_peptides:
            n_samples = len(normalize_df.columns)
            normalize_df = normalize_df.dropna(thresh=round(n_samples * 0.5))

        # Data normalization using the model defined by the user.
        normalized_matrix = scaler.fit_transform(normalize_df)

        # Imputation of the values using KNNImputer (https://scikit-learn.org/0.16/modules/generated/sklearn.preprocessing.Imputer.html)
        if imputation:
           imputer = KNNImputer(n_neighbors=2, weights="uniform")
           normalized_matrix = imputer.fit_transform(normalize_df)

        normalize_df[:] = normalized_matrix
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(id_vars=[PEPTIDE_CHARGE, PEPTIDE_SEQUENCE, CONDITION], var_name=[INTENSITY,class_field], value_vars=INTENSITY)
        normalize_df.pop(INTENSITY)
        normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
        dataset = pd.merge(dataset, normalize_df, how='left', on=[PEPTIDE_CHARGE, PEPTIDE_SEQUENCE, CONDITION, class_field])

    return dataset


@click.command()
@click.option("-p", "--peptides", help="Peptides files from the peptide file generation tool")
@click.option("-c", "--contaminants", help="Contaminants and high abundant proteins to be removed")
@click.option("-l", "--routliers", help="Remove outliers from the peptide table", is_flag=True)
@click.option("-o", "--output", help="Peptide intensity file including other all properties for normalization")
@click.option("-v", "--verbose", help="Print addition information about the distributions of the intensities, number of peptides remove after normalization, etc.", is_flag=True)
def peptide_normalization(peptides: str, contaminants: str, routliers: bool, output: str, verbose: bool) -> None:

    if peptides is None or output is None:
        print_help_msg(peptide_normalization)
        exit(1)

    dataset_df = pd.read_csv(peptides, sep="\t")
    print_dataset_size(dataset_df, "Number of peptides: ", verbose)

    # Print the distribution of the original peptide intensities from quantms analysis
    if verbose:
        plot_distributions(dataset_df, INTENSITY, SAMPLE_ID)
        plot_box_plot(dataset_df, INTENSITY, SAMPLE_ID, log2=True)

    # Remove high abundant and contaminants proteins and the outliers
    if contaminants is not None:
        dataset_df = remove_contaminants_decoys(dataset_df, "contaminants_ids.tsv")
    print_dataset_size(dataset_df, "Peptides after contaminants removal: ", verbose)
    if routliers:
        dataset_df = remove_outliers_iqr(dataset_df)
    print_dataset_size(dataset_df, "Peptides after outliers removal: ", verbose)

    if verbose:
        plot_distributions(dataset_df, INTENSITY, SAMPLE_ID)
        plot_box_plot(dataset_df, INTENSITY, SAMPLE_ID, log2=True)

    dataset_df = intensity_normalization(dataset_df, field=INTENSITY, class_field=SAMPLE_ID, imputation=True, remove_peptides=True, method="robusts")
    if verbose:
        plot_distributions(dataset_df, INTENSITY, SAMPLE_ID)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=True)



if __name__ == '__main__':
    peptide_normalization()