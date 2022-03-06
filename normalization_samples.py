from pandas import DataFrame
from sklearn.preprocessing import MinMaxScaler, StandardScaler, MaxAbsScaler, QuantileTransformer
from sklearn.preprocessing import RobustScaler
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def remove_outliers(dataset: DataFrame):
  """
  Outliers removal
  :param dataset:
  :return:
  """
  Q1 = dataset['Intensity'].quantile(0.25)
  Q3 = dataset['Intensity'].quantile(0.75)
  IQR = Q3 - Q1

  # Filtering Values between Q1-1.5IQR and Q3+1.5IQR
  dataset = dataset.query('(@Q1 - 1.5 * @IQR) <= Intensity <= (@Q3 + 1.5 * @IQR)')
  return dataset


def plot_quantification_box_plot(dataset, method = None, log2 = True, weigth =10, rotation = 45):
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

  plt.figure(figsize=(weigth , 10))
  chart = sns.boxplot(x="SampleID", y="logE", data=normalized, palette="Set2")
  chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation)
  plt.show()
  return dataset

dataset = pd.read_csv("data/PXD008934-Peptide-Intensities.tsv", sep="\t")

contaminants_reader = open("contaminants_ids.tsv", 'r')
contaminants = contaminants_reader.read().split("\n")
contaminants = [cont for cont in contaminants if cont.strip()]
contaminants.append('CONTAMINANTS')

for contaminant in contaminants:
  dataset.drop(index=dataset[dataset['ProteinName'].str.contains(contaminant)].index, inplace=True)

IQR = stats.iqr(dataset['Intensity'], interpolation = 'midpoint')
print(IQR)

dataset = dataset[['ProteinName', 'PeptideSequence', 'Intensity', 'SampleID']]

plot_quantification_box_plot(dataset, weigth = 10, method=None, log2 = True)

dataset = remove_outliers(dataset)
IQR = stats.iqr(dataset['Intensity'], interpolation = 'midpoint')
print(IQR)

plot_quantification_box_plot(dataset, weigth = 10, method=None, log2 = True)
plot_quantification_box_plot(dataset, weigth = 10, method="quantile", log2 = False)
