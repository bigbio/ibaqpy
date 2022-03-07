import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas import DataFrame
from scipy import stats
from sklearn.preprocessing import MinMaxScaler, StandardScaler, MaxAbsScaler, QuantileTransformer
from sklearn.preprocessing import RobustScaler

from ibaqpy_commons import remove_contaminants_decoys, INTENSITY


def remove_outliers(dataset: DataFrame) -> DataFrame:
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


def plot_quantification_box_plot(dataset: DataFrame, method: str = None, log2: bool = True, weigth: int = 10,
                                 rotation: int = 45) -> DataFrame:
    """
    Normalize using different sklearn normalizers
    :param dataset: dataframe
    :param method: method to be used
    :param log2: scale or not the results
    :param weigth: size of plot
    :param rotation: rotation of the plot
    :return:
    """
    normalized = dataset.copy()

    if method is None:
        normalized['normalized'] = normalized['Intensity']
    else:
        scaler = QuantileTransformer(output_distribution="normal")
        if method is "robusts":
            scaler = RobustScaler()
        if method is "minmax":
            scaler = MinMaxScaler()
        if method is "standard":
            scaler = StandardScaler()
        if method is 'maxabs':
            scaler = MaxAbsScaler()
        normalized['normalized'] = scaler.fit_transform(normalized[['Intensity']])

    np.seterr(divide='ignore')
    if log2:
        normalized['logE'] = np.log2(normalized['normalized'])
    else:
        normalized['logE'] = normalized['normalized']

    plt.figure(figsize=(weigth, 10))
    chart = sns.boxplot(x="SampleID", y="logE", data=normalized, palette="Set2")
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation)
    plt.show()
    return dataset


dataset = pd.read_csv("data/PXD008934-Peptide-Intensities.tsv", sep="\t")

dataset = remove_contaminants_decoys(dataset, "contaminants_ids.tsv")
IQR = stats.iqr(dataset['Intensity'], interpolation='midpoint')
print(IQR)

dataset = dataset[['ProteinName', 'PeptideSequence', 'Intensity', 'SampleID']]

plot_quantification_box_plot(dataset, weigth=10, method=None, log2=True)

dataset = remove_outliers(dataset)
IQR = stats.iqr(dataset['Intensity'], interpolation='midpoint')
print(IQR)

plot_quantification_box_plot(dataset, weigth=10, method=None, log2=True)
plot_quantification_box_plot(dataset, weigth=10, method="quantile", log2=False)
