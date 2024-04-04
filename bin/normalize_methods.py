import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform

def normalize_run(df,sdrf_path,method):
    reps = get_replicate(sdrf_path)
    if(reps>1):
        samples = df['SampleID'].unique()
        for sample in samples:
            runs = df.loc[df['SampleID']==sample,'Run'].unique().tolist()
            if(len(runs)>1):
                sample_df = df.loc[df['SampleID']==sample,:]
                map,base = get_normalize_args(sample_df,runs,method)
                for run in runs:
                    run = str(run)
                    run_intensity = df.loc[(df['SampleID']==sample)&(df['Run']==run),'NormIntensity']
                    df.loc[(df['SampleID']==sample)&(df['Run']==run),'NormIntensity'] = run_intensity / (map[run] / base)
        return df
    else:
        return df
        
def get_replicate(sdrf_path):
    sdrf = pd.read_csv(sdrf_path,sep="\t")
    reps = len(sdrf["comment[technical replicate]"].unique())
    return reps

def get_normalize_args(df,runs,method):
    match method:
        case 'mean':
            return normalize_mean(df,runs)
        case 'median':
            return normalize_median(df,runs)
        case 'iqr':
            return normalize_q(df,runs)
            
def normalize_mean(df,runs):
    map = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = df.loc[df['Run']==run,'NormIntensity'].mean()
        map[run] = run_m
        total += run_m
    avg = total / len(runs)
    return map,avg

def normalize_median(df,runs):
    map = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = df.loc[df['Run']==run,'NormIntensity'].median()
        map[run] = run_m
        total += run_m
    med = total / len(runs)
    return map,med

def normalize_q(df,runs):
    map = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = df.loc[df['Run']==run,'NormIntensity'].quantile([0.75,0.25],interpolation='linear').mean()
        map[run] = run_m
        total += run_m
    q = total / len(runs)
    return map,q

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
        case _:
            return -1

# mean
def mean_normalize(df):
    return df / df.mean()

# median
def median_normalize(df):
    return df / df.median()

#max
def max_normalize(df):
    return df / df.max()

#global
def global_normalize(df):
    return df / df.sum()

#max-min
def max_min_mormalize(df):
    min = df.min()
    return (df - min) / (df.max() - min)

#quantile
def quantile_normalize(df):
    index = df.index
    columns = df.columns
    DF = quantile_transform(df)
    df = pd.DataFrame(df,columns=columns,index=index)
    return df