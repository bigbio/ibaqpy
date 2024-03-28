import numpy as np
import pandas as pd
from sklearn.preprocessing import robust_scale
from sklearn.preprocessing import power_transform
from sklearn.preprocessing import quantile_transform

def normalize(df,method):
    match method:
        case 'mean':
            return mean_normalize(df)
        case 'median':
            return median_normalize(df)
        case 'max':
            return max_normalize(df)
        case 'global':
            return global_normalize(df)
        case 'max_min':
            return max_min_mormalize(df)
        case 'z_score':
            return z_score_normalize(df)
        case 'iqr':
            return iqr_normalize(df)
        case 'robust':
            return robust_normalize(df)
        case 'vsn':
            return vsn_normalize(df)
        case 'quantile':
            return quantile_normalize(df)
        case _:
            return -1

# mean
def mean_normalize(df):
    return df / df.mean(axis=0)

# median
def median_normalize(df):
    return df / df.median(axis=0)

#max
def max_normalize(df):
    return df / df.max(axis=0)

#global
def global_normalize(df):
    return df / df.sum(axis=0)

#max-min
def max_min_mormalize(df):
    min = df.min(axis=0)
    return (df - min) / (df.max(axis=0) - min)

#z-score 
def z_score_normalize(df):
    return (df - df.mean(axis=0)) / df.var(axis=0)

#IQR
def iqr_normalize(df):
    Q = df.quantile([0.75,0.25],interpolation='linear',axis=0)
    IQR = Q.loc[0.75,:] - Q.loc[0.25,:]
    return (df - df.median(axis=0)) / IQR

#rubust
def robust_normalize(df):
    index = df.index
    columns = df.columns
    df = robust_scale(df, axis=0)
    df = pd.DataFrame(df,columns=columns,index=index)
    return df

#vsn
def vsn_normalize(df):
    index = df.index
    columns = df.columns
    df = power_transform(df, method='box-cox')
    df = pd.DataFrame(df,columns=columns,index=index)
    return df

#quantile
def quantile_normalize(df):
    index = df.index
    columns = df.columns
    DF = quantile_transform(df,axis=0)
    df = pd.DataFrame(df,columns=columns,index=index)
    return df