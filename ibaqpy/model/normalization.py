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
    def from_str(cls, name: str) -> "FeatureNormalizationMethod":
        """
        Get the normalization method from a string.
        Parameters
        ----------
        name: str The name of the normalization method.

        Returns
        -------
        FeatureNormalizationMethod: The normalization method.
        """
        if name is None:
            return cls.NONE
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def register_replicate_fn(
        self, fn: Callable[[pd.Series], pd.Series]
    ) -> Callable[[pd.Series], pd.Series]:
        _method_registry[self] = fn
        return fn

    def normalize_replicates(self, df: pd.DataFrame, *args, **kwargs):
        """
        Normalize the replicate intensities in the given DataFrame using a registered
        normalization function.

        Parameters:
            df (pd.DataFrame): The DataFrame containing replicate intensity data.
            *args: Additional positional arguments for the normalization function.
            **kwargs: Additional keyword arguments for the normalization function.

        Returns:
            pd.Series: The normalized replicate intensities.
        """
        fn = _method_registry[self]
        return fn(df, *args, **kwargs)

    def normalize_sample(self, df, runs: list[str]) -> tuple[dict[str, pd.Series], float]:
        """
        Normalize replicate intensities for a given sample across multiple runs.

        Parameters:
            df (pd.DataFrame): The DataFrame containing replicate intensity data.
            runs (list[str]): A list of run identifiers for the sample.

        Returns:
            tuple[dict[str, pd.Series], float]: A dictionary mapping each run to its
            normalized replicate intensities and the average metric across all runs.
        """
        map_ = {}
        total = 0
        for run in runs:
            run = str(run)
            run_m = self.normalize_replicates(df.loc[df[TECHREPLICATE] == run, NORM_INTENSITY])
            map_[run] = run_m
            total += run_m
        sample_average_metric = total / len(runs)
        return map_, sample_average_metric

    def normalize_runs(self, df: pd.DataFrame, technical_replicates: int):
        """
        Normalize the intensities of runs in the given DataFrame using a registered

        Parameters
        ----------
        df: pd.DataFrame The DataFrame containing replicate intensity data.
        technical_replicates: int The number of technical replicates for each sample.

        Returns
        -------
        pd.DataFrame: The DataFrame with normalized replicate intensities.

        """
        if technical_replicates > 1:
            samples = df[SAMPLE_ID].unique()
            for sample in samples:
                runs = df.loc[df[SAMPLE_ID] == sample, TECHREPLICATE].unique().tolist()
                if len(runs) > 1:
                    sample_df = df.loc[df[SAMPLE_ID] == sample, :]

                    replicate_metric_map, sample_average_metric = self.normalize_sample(
                        sample_df, runs
                    )

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
    """
    No normalization is performed on the data.
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the replicate intensity data.

    """
    return df


@FeatureNormalizationMethod.Mean.register_replicate_fn
def mean_normalize(df, *args, **kwargs):
    """
    Mean normalization of the data.

    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.

    """
    return df / df.mean()


@FeatureNormalizationMethod.Median.register_replicate_fn
def median_normalize(df, *args, **kwargs):
    """
    Median normalization of the data.
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.

    """
    return df / df.median()


@FeatureNormalizationMethod.Max.register_replicate_fn
def max_normalize(df, *args, **kwargs):
    """
    Max normalization of the data.
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.
    """
    return df / df.max()


@FeatureNormalizationMethod.Global.register_replicate_fn
def global_normalize(df, *args, **kwargs):
    """
    Global normalization of the data.
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.
    """
    return df / df.sum()


@FeatureNormalizationMethod.Max_Min.register_replicate_fn
def max_min_normalize(df, *args, **kwargs):
    """
    Max-Min normalization of the data
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.
    """
    min_ = df.min()
    return (df - min_) / (df.max() - min_)


@FeatureNormalizationMethod.IQR.register_replicate_fn
def iqr_normalization(df, *args, **kwargs):
    """
    IQR normalization of the data.
    Parameters:
    df: pd.DataFrame The DataFrame containing replicate intensity data.
    args: Additional positional arguments
    kwargs: Additional keyword arguments

    Returns:
    pd.DataFrame: The DataFrame containing the normalized replicate intensity data.
    """
    return df.quantile([0.75, 0.25], interpolation="linear").mean()


_peptide_method_registry = {}


class PeptideNormalizationMethod(Enum):
    """
    Enum class for peptide normalization methods, providing functionality to register
    and apply normalization functions to peptide data.

    Attributes:
        NONE: No normalization.
        GlobalMedian: Normalization using global median.
        ConditionMedian: Normalization using condition-specific median.

    Methods:
        from_str(name): Converts a string to a PeptideNormalizationMethod.
        register_replicate_fn(fn): Registers a function for a specific normalization method.
        normalize_sample(dataset_df, sample, med_map): Applies the registered normalization
            function to a sample.
        __call__(dataset_df, sample, med_map): Invokes normalize_sample method.
    """

    NONE = auto()

    GlobalMedian = auto()
    ConditionMedian = auto()

    @classmethod
    def from_str(cls, name: str) -> "PeptideNormalizationMethod":
        """
        Converts a string to a PeptideNormalizationMethod.
        Parameters
        ----------
        name: str The name of the normalization method.

        Returns
        -------
        PeptideNormalizationMethod: The normalization method.
        """
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    def register_replicate_fn(
        self, fn: Callable[[pd.DataFrame, str, dict], pd.DataFrame]
    ) -> Callable[[pd.DataFrame, str, dict], pd.DataFrame]:
        """
        Registers a function for a specific normalization method.
        Parameters
        ----------
        fn: Callable[[pd.DataFrame, str, dict], pd.DataFrame] The normalization function.

        Returns
        -------
        Callable[[pd.DataFrame, str, dict], pd.DataFrame]: The normalization function.
        """
        _peptide_method_registry[self] = fn
        return fn

    def normalize_sample(self, dataset_df: pd.DataFrame, sample: str, med_map: dict):
        """
        Applies the registered normalization function to a sample.
        Parameters
        ----------
        dataset_df: pd.DataFrame The DataFrame containing peptide intensity data.
        sample: str The sample identifier.
        med_map: dict The median map.

        Returns
        -------
        pd.DataFrame: The DataFrame containing the normalized peptide intensity data.
        """
        fn = _peptide_method_registry[self]
        return fn(dataset_df, sample, med_map)

    def __call__(self, dataset_df: pd.DataFrame, sample: str, med_map: dict):
        """
        Invokes the normalize_sample method.
        Parameters
        ----------
        dataset_df: pd.DataFrame The DataFrame containing peptide intensity data.
        sample: str The sample identifier.
        med_map: dict The median map.

        Returns
        -------
        pd.DataFrame: The DataFrame containing the normalized peptide intensity data.
        """
        return self.normalize_sample(dataset_df, sample, med_map)


@PeptideNormalizationMethod.GlobalMedian.register_replicate_fn
def global_median(dataset_df, sample: str, med_map: dict):
    """
    Global median normalization of the data.
    Parameters:
    dataset_df: pd.DataFrame The DataFrame containing peptide intensity data.
    sample: str The sample identifier.
    med_map: dict The median map.

    Returns:
    pd.DataFrame: The DataFrame containing the normalized peptide intensity data.
    """
    dataset_df.loc[:, NORM_INTENSITY] = dataset_df[NORM_INTENSITY] / med_map[sample]
    return dataset_df


@PeptideNormalizationMethod.ConditionMedian.register_replicate_fn
def condition_median(dataset_df, sample: str, med_map: dict):
    """
    Condition median normalization of the data.
    Parameters:
    dataset_df: pd.DataFrame The DataFrame containing peptide intensity data.
    sample: str The sample identifier.
    med_map: dict The median map.

    Returns:
    pd.DataFrame: The DataFrame containing the normalized peptide intensity data.
    """
    con = dataset_df[CONDITION].unique()[0]
    dataset_df.loc[:, NORM_INTENSITY] = dataset_df[NORM_INTENSITY] / med_map[con][sample]


@PeptideNormalizationMethod.NONE.register_replicate_fn
def peptide_no_normalization(dataset_df, sample, med_map):
    """
    No normalization is performed on the data.
    Parameters:
    dataset_df: pd.DataFrame The DataFrame containing peptide intensity data.
    sample: str The sample identifier.
    med_map: dict The median map.

    Returns:
    pd.DataFrame: The DataFrame containing the peptide intensity data.
    """
    return dataset_df
