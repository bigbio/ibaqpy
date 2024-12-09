import pandas as pd
from sklearn.preprocessing import quantile_transform
from ibaqpy.ibaq.ibaqpy_commons import SAMPLE_ID, NORM_INTENSITY, TECHREPLICATE


def normalize_run(df, reps, method):
    if reps > 1:
        samples = df[SAMPLE_ID].unique()
        for sample in samples:
            runs = df.loc[df[SAMPLE_ID] == sample, TECHREPLICATE].unique().tolist()
            if len(runs) > 1:
                sample_df = df.loc[df[SAMPLE_ID] == sample, :]
                map_, base = get_normalize_args(sample_df, runs, method)
                for run in runs:
                    run = str(run)
                    run_intensity = df.loc[
                        (df[SAMPLE_ID] == sample) & (df[TECHREPLICATE] == run),
                        NORM_INTENSITY,
                    ]
                    df.loc[
                        (df[SAMPLE_ID] == sample) & (df[TECHREPLICATE] == run),
                        NORM_INTENSITY,
                    ] = run_intensity / (map_[run] / base)
        return df
    else:
        return df


def get_normalize_args(df, runs, method):
    if method == "mean":
        return normalize_mean(df, runs)
    elif method == "median":
        return normalize_median(df, runs)
    elif method == "iqr":
        return normalize_q(df, runs)
    else:
        exit(f"Method {method} not supported!")


def normalize_mean(df, runs):
    map_ = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = df.loc[df[TECHREPLICATE] == run, NORM_INTENSITY].mean()
        map_[run] = run_m
        total += run_m
    avg = total / len(runs)
    return map_, avg


def normalize_median(df, runs):
    map_ = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = df.loc[df[TECHREPLICATE] == run, NORM_INTENSITY].median()
        map_[run] = run_m
        total += run_m
    med = total / len(runs)
    return map_, med


def normalize_q(df, runs):
    map_ = {}
    total = 0
    for run in runs:
        run = str(run)
        run_m = (
            df.loc[df[TECHREPLICATE] == run, NORM_INTENSITY]
            .quantile([0.75, 0.25], interpolation="linear")
            .mean()
        )
        map_[run] = run_m
        total += run_m
    q = total / len(runs)
    return map_, q


def normalize(df, method):
    if method == "mean":
        return mean_normalize(df)
    elif method == "median":
        return median_normalize(df)
    elif method == "max":
        return max_normalize(df)
    elif method == "global":
        return global_normalize(df)
    elif method == "max_min":
        return max_min_mormalize(df)
    else:
        exit(f"Method {method} not supported!")


# mean
def mean_normalize(df):
    return df / df.mean()


# median
def median_normalize(df):
    return df / df.median()


# max
def max_normalize(df):
    return df / df.max()


# global
def global_normalize(df):
    return df / df.sum()


# max-min
def max_min_mormalize(df):
    min_ = df.min()
    return (df - min_) / (df.max() - min_)


# quantile
def quantile_normalize(df):
    index = df.index
    columns = df.columns
    df = quantile_transform(df)
    df = pd.DataFrame(df, columns=columns, index=index)
    return df
