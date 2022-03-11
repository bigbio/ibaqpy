import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas import DataFrame
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.preprocessing import MinMaxScaler, StandardScaler, MaxAbsScaler, QuantileTransformer, Normalizer
from sklearn.preprocessing import RobustScaler
import qnorm

from ibaqpy_commons import remove_contaminants_decoys, INTENSITY, SAMPLE_ID, NORM_INTENSITY, \
    PEPTIDE_SEQUENCE, CONDITION, PEPTIDE_CHARGE, FRACTION, RUN, BIOREPLICATE, RT


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
    Q1 = dataset[INTENSITY].quantile(0.25)
    Q3 = dataset[INTENSITY].quantile(0.75)
    IQR = Q3 - Q1

    dataset.query('(@Q1 - 1.5 * @IQR) <= Intensity <= (@Q3 + 1.5 * @IQR)', inplace=True)


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


def intensity_normalization(dataset: DataFrame, field: str, class_field: str = "all", scaling_method: str = "msstats",
                            imputation: bool = False,
                            imputation_method: str = "simple", remove_peptides: bool = False) -> DataFrame:

    ## TODO add imputation and/or removal to those two norm strategies
    if scaling_method == 'msstats':
        g = dataset.groupby(['Run', 'Fraction'])[INTENSITY].apply(np.median)
        g.name = 'RunMedian'
        dataset = dataset.join(g, on=['Run', 'Fraction'])
        # TODO might be quicker with transform but I could not make it work in short time for multidim. groupby
        #dataset['RunMedian'] = dataset[INTENSITY].groupby(dataset[['Run', 'Fraction']]).transform('median')
        dataset['FractionMedian'] = dataset['RunMedian'].groupby(dataset['Fraction']).transform('median')
        dataset[NORM_INTENSITY] = dataset[INTENSITY] - dataset['RunMedian'] + dataset['FractionMedian']
        return dataset

    elif scaling_method == 'qnorm':
        # pivot to have one col per sample
        normalize_df = pd.pivot_table(dataset, index=[PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, FRACTION, RUN, BIOREPLICATE, RT],
                                      columns=class_field, values=field, aggfunc={field: np.mean})
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        print(normalize_df.head())
        dataset = normalize_df.melt(value_name=NORM_INTENSITY)
        print(dataset.head())
        return dataset

    # Scaler selection
    elif scaling_method == "standard":
        scaler = StandardScaler(with_mean=False)
    elif scaling_method == "quantile":
        scaler = QuantileTransformer(output_distribution="normal")
    elif scaling_method == "normalizer":
        scaler = Normalizer()
    elif scaling_method == "robust":
        scaler = RobustScaler()
    elif scaling_method == "minmax":
        scaler = MinMaxScaler()
    elif scaling_method == 'maxabs':
        scaler = MaxAbsScaler()

    # normalize the intensities across all samples
    if class_field == "all":
        normalize_df = dataset[[field]]
        dataset[NORM_INTENSITY] = scaler.fit_transform(normalize_df)
    else:
        # normalize taking into account samples
        normalize_df = dataset[[PEPTIDE_SEQUENCE, CONDITION, field, class_field]]
        # group peptide + charge into a single peptide intensity using the mean.
        normalize_df = pd.pivot_table(normalize_df, values=field, index=[PEPTIDE_SEQUENCE, CONDITION],
                                      columns=class_field,
                                      aggfunc={field: np.mean})

        print(normalize_df.head())

        # Remove all peptides in less than 30% of the samples.
        if remove_peptides:
            n_samples = len(normalize_df.columns)
            normalize_df = normalize_df.dropna(thresh=round(n_samples * 0.3))

        # Imputation of the values using KNNImputer
        # (https://scikit-learn.org/0.16/modules/generated/sklearn.preprocessing.Imputer.html)
        if imputation:
            imputer = SimpleImputer()
            if imputation_method != "simple":
                imputer = KNNImputer(n_neighbors=2, weights="uniform")
            normalized_matrix = imputer.fit_transform(normalize_df)
            if scaling_method != "quantile":
                normalized_matrix = scaler.fit_transform(normalized_matrix)
            else:
                normalized_matrix = qnorm.quantile_normalize(normalized_matrix)
        else:
            normalized_matrix = scaler.fit_transform(normalize_df)

        normalize_df[:] = normalized_matrix
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(id_vars=[PEPTIDE_SEQUENCE, CONDITION])
        normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
        dataset = pd.merge(dataset, normalize_df, how='left', on=[PEPTIDE_SEQUENCE, CONDITION, class_field])
        dataset.pop(NORM_INTENSITY + "_x")
        dataset.rename(
            columns={NORM_INTENSITY + "_y": NORM_INTENSITY}, inplace=True)

    return dataset


@click.command()
@click.option("--peptides", help="Peptides files from the peptide file generation tool")
@click.option("--contaminants", help="Contaminants and high abundant proteins to be removed")
@click.option("--routliers", help="Remove outliers from the peptide table", is_flag=True)
@click.option("--output", help="Peptide intensity file including other all properties for normalization")
@click.option('--nmethod', help="Normalization method used to normalize intensities for all samples (options: quantile, robusts, standard)", default="quantile")
@click.option("--compress", help="Read the input peptides file in compress gzip file", is_flag= True)
@click.option("--log2", help="Transform to log2 the peptide intensity values before normalization", is_flag=True)
@click.option("--violin", help="Use violin plot instead of boxplot for distribution representations", is_flag=True)
@click.option("--verbose",
              help="Print addition information about the distributions of the intensities, number of peptides remove after normalization, etc.",
              is_flag=True)
def peptide_normalization(peptides: str, contaminants: str, routliers: bool, output: str, nmethod: str, compress: bool, log2: bool,
                          violin: bool, verbose: bool) -> None:
    if peptides is None or output is None:
        print_help_msg(peptide_normalization)
        exit(1)

    pd.set_option('display.max_columns', None)
    # TODO infer from filename
    print("Loading data..")
    compression_method = 'gzip' if compress else None
    if compress:
        dataset_df = pd.read_csv(peptides, sep="\t", compression=compression_method)
    else:
        dataset_df = pd.read_csv(peptides, sep="\t")
    print_dataset_size(dataset_df, "Number of peptides: ", verbose)

    print("Logarithmize if specified..")
    dataset_df[NORM_INTENSITY] = np.log2(dataset_df[INTENSITY]) if log2 else dataset_df[INTENSITY]

    # Print the distribution of the original peptide intensities from quantms analysis
    if verbose:
        plot_distributions(dataset_df, INTENSITY, SAMPLE_ID, log2=not log2)
        plot_box_plot(dataset_df, INTENSITY, SAMPLE_ID, log2=not log2,
                      title="Original peptide intensity distribution (no normalization)", violin=violin)

    # Remove high abundant and contaminants proteins and the outliers
    if contaminants is not None:
        print("Remove contaminants...")
        dataset_df = remove_contaminants_decoys(dataset_df, "contaminants_ids.tsv")
    print_dataset_size(dataset_df, "Peptides after contaminants removal: ", verbose)

    if verbose:
        plot_distributions(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=not log2)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=not log2,
                      title="Peptide intensity distribution after contaminants removal", violin=violin)

    dataset_df = intensity_normalization(dataset_df, field=NORM_INTENSITY, class_field=SAMPLE_ID,
                                         imputation=True, remove_peptides=routliers, scaling_method=nmethod)
    if verbose:
        log_after_norm = nmethod == "msstats" or nmethod == "qnorm" or ((nmethod == "quantile" or nmethod == "robust") and not log2)
        plot_distributions(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm,
                      title="Peptide intensity distribution after imputation, normalization", violin=violin)


if __name__ == '__main__':
    peptide_normalization()
