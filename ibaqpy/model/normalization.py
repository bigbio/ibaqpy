from typing import Callable
from enum import Enum, auto

import pandas as pd

from ibaqpy.ibaq.ibaqpy_commons import CONDITION, NORM_INTENSITY, SAMPLE_ID, TECHREPLICATE

_method_registry: dict["FeatureNormalizationMethod", Callable[[pd.Series], pd.Series]] = {}

class FeatureNormalizationMethod(Enum):
    NONE = auto()

    Mean = auto()
    Median = auto()
    Max = auto()
    Global = auto()
    Max_Min = auto()
    IQR = auto()

    @classmethod
    def from_str(cls, name: str) -> 'FeatureNormalizationMethod':
        if name is None:
            return cls.NONE
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def register_replicate_fn(self, fn: Callable[[pd.Series], pd.Series]) -> Callable[[pd.Series], pd.Series]:
        _method_registry[self] = fn
        return fn

    def normalize_replicates(self, df: pd.DataFrame, *args, **kwargs):
        fn = _method_registry[self]
        return fn(df, *args, **kwargs)

    def normalize_sample(self, df, runs: list[str]) -> tuple[dict[str, pd.Series], float]:
        map_ = {}
        total = 0
        for run in runs:
            run = str(run)
            run_m = (
                self.normalize_replicates(df.loc[df[TECHREPLICATE] == run, NORM_INTENSITY])
            )
            map_[run] = run_m
            total += run_m
        sample_average_metric = total / len(runs)
        return map_, sample_average_metric

    def normalize_runs(self, df: pd.DataFrame, technical_replicates: int):
        if technical_replicates > 1:
            samples = df[SAMPLE_ID].unique()
            for sample in samples:
                runs = df.loc[df[SAMPLE_ID] == sample, TECHREPLICATE].unique().tolist()
                if len(runs) > 1:
                    sample_df = df.loc[df[SAMPLE_ID] == sample, :]

                    replicate_metric_map, sample_average_metric = self.normalize_sample(sample_df, runs)

                    # For each replicate in each sample, normalize the per-replicate
                    # intensity by a replicate-level statistic, relative to the sample
                    # average over that replicate statistic.
                    #
                    # In effect, this scales runs down when the replicate average > sample average
                    # and scales runs up when the replicate average < sample average.
                    for run in runs:
                        run = str(run)
                        run_intensity = df.loc[
                            (df[SAMPLE_ID] == sample) & (df[TECHREPLICATE] == run),
                            NORM_INTENSITY,
                        ]
                        df.loc[
                            (df[SAMPLE_ID] == sample) & (df[TECHREPLICATE] == run),
                            NORM_INTENSITY,
                        ] = run_intensity / (replicate_metric_map[run] / sample_average_metric)
            return df
        else:
            return df

    def __call__(self, df: pd.DataFrame, technical_replicates: int):
        return self.normalize_runs(df, technical_replicates)

@FeatureNormalizationMethod.NONE.register_replicate_fn
def no_normalization(df, *args, **kwargs):
    return df


@FeatureNormalizationMethod.Mean.register_replicate_fn
def mean_normalize(df, *args, **kwargs):
    return df / df.mean()


@FeatureNormalizationMethod.Median.register_replicate_fn
def median_normalize(df, *args, **kwargs):
    return df / df.median()


@FeatureNormalizationMethod.Max.register_replicate_fn
def max_normalize(df, *args, **kwargs):
    return df / df.max()


@FeatureNormalizationMethod.Global.register_replicate_fn
def global_normalize(df, *args, **kwargs):
    return df / df.sum()


@FeatureNormalizationMethod.Max_Min.register_replicate_fn
def max_min_normalize(df, *args, **kwargs):
    min_ = df.min()
    return (df - min_) / (df.max() - min_)


@FeatureNormalizationMethod.IQR.register_replicate_fn
def iqr_normalization(df, *args, **kwargs):
    return df.quantile([0.75, 0.25], interpolation="linear").mean()


_peptide_method_registry = {}


class PeptideNormalizationMethod(Enum):
    NONE = auto()

    GlobalMedian = auto()
    ConditionMedian = auto()

    @classmethod
    def from_str(cls, name: str) -> "PeptideNormalizationMethod":
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def register_replicate_fn(
        self, fn: Callable[[pd.DataFrame, str, dict], pd.DataFrame]
    ) -> Callable[[pd.DataFrame, str, dict], pd.DataFrame]:
        _peptide_method_registry[self] = fn
        return fn

    def normalize_sample(self, dataset_df: pd.DataFrame, sample: str, med_map: dict):
        fn = _peptide_method_registry[self]
        return fn(dataset_df, sample, med_map)

    def __call__(self, dataset_df: pd.DataFrame, sample: str, med_map: dict):
        return self.normalize_sample(dataset_df, sample, med_map)


@PeptideNormalizationMethod.GlobalMedian.register_replicate_fn
def global_median(dataset_df, sample: str, med_map: dict):
    dataset_df.loc[:, NORM_INTENSITY] = dataset_df[NORM_INTENSITY] / med_map[sample]
    return dataset_df


@PeptideNormalizationMethod.ConditionMedian.register_replicate_fn
def condition_median(dataset_df, sample: str, med_map: dict):
    con = dataset_df[CONDITION].unique()[0]
    dataset_df.loc[:, NORM_INTENSITY] = (
        dataset_df[NORM_INTENSITY] / med_map[con][sample]
    )


@PeptideNormalizationMethod.NONE.register_replicate_fn
def peptide_no_normalization(dataset_df, sample, med_map):
    return dataset_df