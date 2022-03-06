from pandas import DataFrame
from sklearn.preprocessing import MinMaxScaler, StandardScaler, MaxAbsScaler, QuantileTransformer
from sklearn.preprocessing import RobustScaler
import re
import pandas as pd
import seaborn as sns
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def get_sample(abundance: str):
  m = re.search(r"\[([A-Za-z0-9_]+)\]", abundance)
  sample = m.group(1)

  return sample

def remove_outliers(dataset: DataFrame):
  Q1 = np.percentile(dataset['Intensity'], 25, interpolation='midpoint')

  Q3 = np.percentile(dataset['Intensity'], 75, interpolation='midpoint')
  IQR = Q3 - Q1

  print("Old Shape: ", dataset.shape)

  # Upper bound
  upper = np.where(dataset['Intensity'] >= (Q3 + 1.5 * IQR))
  # Lower bound
  lower = np.where(dataset['Intensity'] <= (Q1 - 1.5 * IQR))

  ''' Removing the Outliers '''
  dataset.drop(upper[-1], inplace=True)
  dataset.drop(lower[-1], inplace=True)
  return dataset


def plot_quantification_box_plot(dataset, method = None, log2 = True, weigth =10, rotation = 45):

  normalized = dataset.copy()

  if method is None:
    normalized['normalized'] = normalized['Intensity']
  if method is "robusts":
    scaler = RobustScaler()
    normalized['normalized'] = scaler.fit_transform(normalized[['Intensity']])
  if method is "minmax":
    scaler = MinMaxScaler()
    normalized['normalized'] = scaler.fit_transform(normalized[['Intensity']])
  if method is "standard":
    scaler = StandardScaler()
    normalized['normalized'] = scaler.fit_transform(normalized[['Intensity']])

  if method is 'maxabs':
    scaler = MaxAbsScaler()
    normalized['normalized'] = scaler.fit_transform(normalized[['Intensity']])

  if method is 'quantile':
    scaler = QuantileTransformer(output_distribution="normal")
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
