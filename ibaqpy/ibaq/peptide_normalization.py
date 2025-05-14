import os
import re
import time
from typing import Iterator, Optional

import pandas as pd
import numpy as np
import duckdb

from ibaqpy.model.quantification_type import QuantificationCategory, IsobaricLabel
from ibaqpy.model.normalization import FeatureNormalizationMethod, PeptideNormalizationMethod
from ibaqpy.ibaq.ibaqpy_commons import (
    BIOREPLICATE,
    TECHREPLICATE,
    CHANNEL,
    CONDITION,
    FRACTION,
    INTENSITY,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PEPTIDE_CHARGE,
    PEPTIDE_SEQUENCE,
    PROTEIN_NAME,
    RUN,
    SAMPLE_ID,
    PARQUET_COLUMNS,
    parquet_map,
)

from .write_queue import WriteParquetTask, WriteCSVTask
from .logger import get_logger, log_execution_time

# Get a logger for this module
logger = get_logger("ibaqpy.peptide_normalization")


def parse_uniprot_accession(uniprot_id: str) -> str:
    """
    Parse a UniProt accession string to extract and return the core accession numbers.

    This function takes a UniProt ID string, which may contain multiple accessions
    separated by semicolons. Each accession may have a format with two pipe ('|')
    characters, in which case the core accession number is extracted. The function
    returns a semicolon-separated string of the core accession numbers.

    Parameters:
        uniprot_id (str): A string containing one or more UniProt accessions.

    Returns:
        str: A semicolon-separated string of core accession numbers.
    """
    uniprot_list = uniprot_id.split(";")
    result_uniprot_list = []
    for accession in uniprot_list:
        if accession.count("|") == 2:
            accession = accession.split("|")[1]
        result_uniprot_list.append(accession)
    return ";".join(result_uniprot_list)


def get_canonical_peptide(peptide_sequence: str) -> str:
    """
    Remove modifications and special characters from a peptide sequence.

    This function takes a peptide sequence string and removes any modifications
    enclosed in parentheses or brackets, as well as any periods or hyphens,
    returning the cleaned canonical peptide sequence.

    Parameters:
        peptide_sequence (str): The peptide sequence to be cleaned.

    Returns:
        str: The cleaned canonical peptide sequence.
    """
    clean_peptide = re.sub(r"[\(\[].*?[\)\]]", "", peptide_sequence)
    clean_peptide = clean_peptide.replace(".", "").replace("-", "")
    return clean_peptide


def analyse_sdrf(
    sdrf_path: str,
) -> tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
    """
    Analyzes an SDRF file to determine quantification details.

    Parameters:
        sdrf_path (str): The file path to the SDRF file.

    Returns:
        tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
        A tuple containing the number of technical repetitions, the quantification category,
        a list of unique sample names, and the isobaric label scheme if applicable.
    """
    sdrf_df = pd.read_csv(sdrf_path, sep="\t")
    sdrf_df.columns = [i.lower() for i in sdrf_df.columns]

    labels = set(sdrf_df["comment[label]"])
    # Determine label type
    label, channel_set = QuantificationCategory.classify(labels)
    if label in (QuantificationCategory.TMT, QuantificationCategory.ITRAQ):
        choice_df = (
            pd.DataFrame.from_dict(channel_set.channels(), orient="index", columns=[CHANNEL])
            .reset_index()
            .rename(columns={"index": "comment[label]"})
        )
        sdrf_df = sdrf_df.merge(choice_df, on="comment[label]", how="left")
    sample_names = sdrf_df["source name"].unique().tolist()
    technical_repetitions = len(sdrf_df["comment[technical replicate]"].unique())
    return technical_repetitions, label, sample_names, channel_set


def remove_contaminants_entrapments_decoys(
    dataset: pd.DataFrame, protein_field=PROTEIN_NAME
) -> pd.DataFrame:
    """
    Remove rows from the dataset that contain contaminants, entrapments, or decoys.

    This function filters out entries in the specified protein field that match
    any of the terms "CONTAMINANT", "ENTRAP", or "DECOY".

    Parameters:
        dataset (pd.DataFrame): The input DataFrame containing protein data.
        protein_field (str): The column name in the DataFrame to check for contaminants.

    Returns:
        pd.DataFrame: A DataFrame with the contaminants, entrapments, and decoys removed.
    """
    contaminants = ["CONTAMINANT", "ENTRAP", "DECOY"]
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex)]


def remove_protein_by_ids(
    dataset: pd.DataFrame, protein_file: str, protein_field=PROTEIN_NAME
) -> pd.DataFrame:
    """
    Remove proteins from a dataset based on a list of protein IDs.

    This function reads a file containing protein IDs and removes any rows
    from the dataset where the specified protein field matches any of the IDs.

    Parameters:
        dataset (pd.DataFrame): The dataset containing protein information.
        protein_file (str): Path to the file containing protein IDs to be removed.
        protein_field (str): The field in the dataset to check for protein IDs.
                             Defaults to PROTEIN_NAME.

    Returns:
        pd.DataFrame: A DataFrame with the specified proteins removed.
    """
    contaminants_reader = open(protein_file, "r")
    contaminants = contaminants_reader.read().split("\n")
    contaminants = [cont for cont in contaminants if cont.strip()]
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex, regex=True)]


def reformat_quantms_feature_table_quant_labels(
    data_df: pd.DataFrame, label: QuantificationCategory, choice: Optional[IsobaricLabel]
) -> pd.DataFrame:
    """
    Reformats a DataFrame containing quantification labels for QuantMS features.

    This function renames columns in the input DataFrame according to a predefined mapping
    and processes the protein names and channel information based on the quantification
    category and isobaric label choice.

    Parameters:
        data_df (pd.DataFrame): The input DataFrame containing quantification data.
        label (QuantificationCategory): The quantification category (e.g., LFQ, TMT, ITRAQ).
        choice (Optional[IsobaricLabel]): The isobaric label scheme, if applicable.

    Returns:
        pd.DataFrame: The reformatted DataFrame with updated column names and channel information.
    """
    data_df = data_df.rename(columns=parquet_map)
    data_df[PROTEIN_NAME] = data_df[PROTEIN_NAME].str.join(";")
    if label == QuantificationCategory.LFQ:
        data_df.drop(CHANNEL, inplace=True, axis=1)
    else:
        data_df[CHANNEL] = data_df[CHANNEL].map(choice.channels())

    return data_df


def apply_initial_filtering(data_df: pd.DataFrame, min_aa: int) -> pd.DataFrame:
    """
    Apply initial filtering to a DataFrame containing peptide data.

    This function filters out rows with zero intensity, removes entries with
    'Empty' conditions unless they are null, and excludes peptides with fewer
    amino acids than the specified minimum. It also processes protein names
    to extract core accession numbers and ensures the presence of a 'Fraction'
    column. Additionally, it extracts technical replicate information from the
    'Run' column if applicable, and retains only relevant columns for further
    analysis.

    Parameters:
        data_df (pd.DataFrame): The input DataFrame containing peptide data.
        min_aa (int): The minimum number of amino acids required for peptides.

    Returns:
        pd.DataFrame: The filtered DataFrame with relevant columns.
    """

    # Remove 0 intensity signals from the data
    data_df = data_df[data_df[INTENSITY] > 0]

    data_df = data_df[(data_df["Condition"] != "Empty") | (data_df["Condition"].isnull())]

    # Filter peptides with less amino acids than min_aa (default: 7)
    data_df.loc[:, "len"] = data_df[PEPTIDE_CANONICAL].apply(len)
    data_df = data_df[data_df["len"] >= min_aa]
    data_df.drop(["len"], inplace=True, axis=1)
    data_df[PROTEIN_NAME] = data_df[PROTEIN_NAME].apply(parse_uniprot_accession)
    if FRACTION not in data_df.columns:
        data_df[FRACTION] = 1

    # TODO: What if there's a mix of identifiers?
    if data_df[RUN].str.contains("_").all():
        data_df[TECHREPLICATE] = data_df[RUN].str.split("_").str.get(1)
        data_df[TECHREPLICATE] = data_df[TECHREPLICATE].astype("int")
    else:
        data_df[TECHREPLICATE] = data_df[RUN].astype("int")

    data_df = data_df[
        [
            PROTEIN_NAME,
            PEPTIDE_SEQUENCE,
            PEPTIDE_CANONICAL,
            PEPTIDE_CHARGE,
            INTENSITY,
            CONDITION,
            TECHREPLICATE,
            BIOREPLICATE,
            FRACTION,
            SAMPLE_ID,
        ]
    ]
    data_df[CONDITION] = pd.Categorical(data_df[CONDITION])
    data_df[SAMPLE_ID] = pd.Categorical(data_df[SAMPLE_ID])

    return data_df


def merge_fractions(dataset: pd.DataFrame) -> pd.DataFrame:
    """
    Merge fractions in the dataset by grouping and aggregating normalized intensity.

    This function removes rows with missing normalized intensity values and groups
    the dataset by protein name, peptide sequence, canonical peptide, peptide charge,
    condition, biological replicate, technical replicate, and sample ID. It then
    aggregates the normalized intensity by taking the maximum value for each group.

    Parameters:
        dataset (pd.DataFrame): The input DataFrame containing peptide data.

    Returns:
        pd.DataFrame: A DataFrame with merged fractions and maximum normalized intensity.
    """
    dataset.dropna(subset=[NORM_INTENSITY], inplace=True)
    dataset = dataset.groupby(
        [
            PROTEIN_NAME,
            PEPTIDE_SEQUENCE,
            PEPTIDE_CANONICAL,
            PEPTIDE_CHARGE,
            CONDITION,
            BIOREPLICATE,
            TECHREPLICATE,
            SAMPLE_ID,
        ],
        observed=True,
    ).agg({NORM_INTENSITY: "max"})
    dataset.reset_index(inplace=True)
    return dataset


def get_peptidoform_normalize_intensities(
    dataset: pd.DataFrame, higher_intensity: bool = True
) -> pd.DataFrame:
    """
    Normalize peptide intensities in a dataset by selecting the highest intensity
    for each unique combination of peptide sequence, charge, sample ID, condition,
    and biological replicate.

    Parameters:
        dataset (pd.DataFrame): The input DataFrame containing peptide data.
        higher_intensity (bool): If True, selects the row with the highest normalized
                                 intensity for each group. Defaults to True.

    Returns:
        pd.DataFrame: A DataFrame with normalized intensities, filtered to retain
                      only the highest intensity entries per group.
    """
    dataset.dropna(subset=[NORM_INTENSITY], inplace=True)
    if higher_intensity:
        dataset = dataset.loc[
            dataset.groupby(
                [PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE],
                observed=True,
            )[NORM_INTENSITY].idxmax()
        ]
    # else:
    #     dataset = dataset.loc[
    #         dataset.groupby(
    #             [PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE],
    #             observed=True,
    #         )[SEARCH_ENGINE].idxmax()
    #     ]
    dataset.reset_index(drop=True, inplace=True)
    return dataset


def sum_peptidoform_intensities(dataset: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate normalized intensities for each unique peptidoform.

    This function processes a dataset by summing the normalized intensities
    for each unique combination of protein name, peptidoform, sample ID,
    biological replicate, and condition. It removes any rows with missing
    normalized intensity values, performs the aggregation, and returns a
    DataFrame with unique entries.

    Parameters:
        dataset (pd.DataFrame): The input DataFrame containing peptidoform
        data with normalized intensities.

    Returns:
        pd.DataFrame: A DataFrame with summed normalized intensities for each
        unique peptidoform entry.
    """
    dataset.dropna(subset=[NORM_INTENSITY], inplace=True)
    dataset = dataset[
        [
            PROTEIN_NAME,
            PEPTIDE_CANONICAL,
            SAMPLE_ID,
            BIOREPLICATE,
            CONDITION,
            NORM_INTENSITY,
        ]
    ]
    dataset.loc[:, NORM_INTENSITY] = dataset.groupby(
        [PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION],
        observed=True,
    )[NORM_INTENSITY].transform("sum")
    dataset = dataset.drop_duplicates()
    dataset.reset_index(inplace=True, drop=True)
    return dataset


class Feature:
    """
    Represents a feature in a proteomics dataset, providing methods for data manipulation
    and analysis using a DuckDB database connection to a Parquet file.

    Attributes:
        labels (Optional[list[str]]): A list of unique labels from the dataset.
        label (Optional[QuantificationCategory]): The quantification category inferred from labels.
        choice (Optional[IsobaricLabel]): The isobaric label scheme inferred from labels.
        technical_repetitions (Optional[int]): The number of unique technical repetitions.

    Methods:
        __init__(database_path: str): Initializes the Feature object and connects to the database.
        standardize_df(df: pd.DataFrame) -> pd.DataFrame: Standardizes column names in a DataFrame.
        experimental_inference() -> tuple: Infers experimental details from the dataset.
        low_frequency_peptides(percentage=0.2) -> tuple: Identifies low-frequency peptides.
        csv2parquet(csv): Converts a CSV file to a Parquet file.
        get_report_from_database(samples: list, columns: list = None) -> pd.DataFrame: Retrieves a report from the database.
        iter_samples(sample_num: int = 20, columns: list = None) -> Iterator: Iterates over samples in batches.
        get_unique_samples() -> list[str]: Retrieves unique sample accessions.
        get_unique_labels() -> list[str]: Retrieves unique channel labels.
        get_unique_tec_reps() -> list[int]: Retrieves unique technical repetition identifiers.
        get_median_map() -> dict[str, float]: Computes a median intensity map for samples.
        get_report_condition_from_database(cons: list, columns: list = None) -> pd.DataFrame: Retrieves a report based on conditions.
        iter_conditions(conditions: int = 10, columns: list = None) -> Iterator: Iterates over conditions in batches.
        get_unique_conditions() -> list[str]: Retrieves unique experimental conditions.
        get_median_map_to_condition() -> dict[str, dict[str, float]]: Computes a median intensity map for conditions.
    """

    labels: Optional[list[str]]
    label: Optional[QuantificationCategory]
    choice: Optional[IsobaricLabel]
    technical_repetitions: Optional[int]

    def __init__(self, database_path: str):
        if os.path.exists(database_path):
            self.parquet_db = duckdb.connect()
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(database_path)
            )
            self.samples = self.get_unique_samples()
        else:
            raise FileNotFoundError(f"the file {database_path} does not exist.")

    @staticmethod
    def standardize_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardizes column names in the given DataFrame.

        Parameters:
            df (pd.DataFrame): The DataFrame with columns to be standardized.

        Returns:
            pd.DataFrame: A DataFrame with standardized column names.
        """
        return df.rename(
            {"protein_accessions": "pg_accessions", "charge": "precursor_charge"}, axis=1
        )

    @property
    def experimental_inference(
        self,
    ) -> tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
        """
        Infers experimental details from the dataset, including the number of unique
        technical repetitions, the quantification category, the list of samples, and
        the isobaric label scheme.

        Returns:
            tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
            A tuple containing the number of unique technical repetitions, the inferred
            quantification category, the list of samples, and the isobaric label scheme.
        """
        self.labels = self.get_unique_labels()
        self.label, self.choice = QuantificationCategory.classify(self.labels)
        self.technical_repetitions = self.get_unique_tec_reps()
        return len(self.technical_repetitions), self.label, self.samples, self.choice

    @property
    def low_frequency_peptides(self, percentage=0.2) -> tuple:
        """
        Identifies peptides that occur with low frequency across samples.

        Parameters:
            percentage (float): The threshold percentage of samples a peptide must appear in
                                to be considered low frequency. Defaults to 0.2.

        Returns:
            tuple: A tuple of tuples, each containing a protein group accession and a peptide sequence
                   that are identified as low frequency.
        """
        f_table = self.parquet_db.sql(
            """
            SELECT "sequence", "pg_accessions", COUNT(DISTINCT sample_accession) as "count" from parquet_db
            GROUP BY "sequence","pg_accessions"
            """
        ).df()
        f_table.dropna(subset=["pg_accessions"], inplace=True)
        try:
            f_table["pg_accessions"] = f_table["pg_accessions"].apply(lambda x: x[0].split("|")[1])
        except IndexError:
            f_table["pg_accessions"] = f_table["pg_accessions"].apply(lambda x: x[0])
        except Exception as e:
            raise ValueError(
                "Some errors occurred when parsing pg_accessions column in feature parquet!"
            ) from e
        f_table.set_index(["sequence", "pg_accessions"], inplace=True)
        f_table.drop(
            f_table[f_table["count"] >= (percentage * len(self.samples))].index,
            inplace=True,
        )
        f_table.reset_index(inplace=True)
        return tuple(zip(f_table["pg_accessions"], f_table["sequence"]))

    @staticmethod
    def csv2parquet(csv):
        """
        Converts a CSV file to a Parquet file using DuckDB.

        Parameters:
            csv (str): The file path of the CSV file to be converted.

        The converted Parquet file will be saved in the same directory
        with the same name as the CSV file, but with a .parquet extension.
        """
        parquet_path = os.path.splitext(csv)[0] + ".parquet"
        duckdb.read_csv(csv).to_parquet(parquet_path)

    def get_report_from_database(self, samples: list, columns: list = None):
        """
        Retrieves a standardized report from the database for specified samples.

        Parameters:
            samples (list): A list of sample accessions to filter the report.
            columns (list, optional): A list of column names to include in the report.
                                      If None, all columns are included.

        Returns:
            pd.DataFrame: A DataFrame containing the report with standardized column names.
        """
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            """SELECT {} FROM parquet_db WHERE sample_accession IN {}""".format(
                cols, tuple(samples)
            )
        )
        report = database.df()
        return Feature.standardize_df(report)

    def iter_samples(
        self, sample_num: int = 20, columns: list = None
    ) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """
        Iterates over samples in batches, yielding each batch along with its corresponding
        report DataFrame.

        Parameters:
            sample_num (int, optional): The number of samples to include in each batch. Defaults to 20.
            columns (list, optional): A list of column names to include in the report. If None, all columns are included.

        Yields:
            Iterator[tuple[list[str], pd.DataFrame]]: An iterator over tuples, each containing a list of sample accessions
                                                      and a DataFrame with the report for those samples.
        """
        ref_list = [
            self.samples[i : i + sample_num] for i in range(0, len(self.samples), sample_num)
        ]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_samples(self) -> list[str]:
        """
        Retrieves a list of unique sample accessions from the Parquet database.

        Returns:
            list[str]: A list of unique sample accession identifiers.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT sample_accession FROM parquet_db").df()
        return unique["sample_accession"].tolist()

    def get_unique_labels(self) -> list[str]:
        """
        Retrieves a list of unique channel labels from the Parquet database.

        Returns:
            list[str]: A list of unique channel labels.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT channel FROM parquet_db").df()
        return unique["channel"].tolist()

    def get_unique_tec_reps(self) -> list[int]:
        """
        Retrieves a list of unique technical repetition identifiers from the Parquet database.

        Returns:
            list[int]: A list of unique technical repetition identifiers as integers.

        Raises:
            ValueError: If there is an error converting the 'run' identifiers to integers.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT run FROM parquet_db").df()
        try:
            if unique["run"].str.contains("_").all():
                unique["run"] = unique["run"].str.split("_").str.get(1)
                unique["run"] = unique["run"].astype("int")
            else:
                unique["run"] = unique["run"].astype("int")
        except ValueError as e:
            raise ValueError(
                f"Some errors occurred when getting technical repetitions: {e}"
            ) from e

        return unique["run"].tolist()

    def get_median_map(self) -> dict[str, float]:
        """
        Computes a median intensity map for samples, normalizing each sample's median
        intensity by the global median intensity across all samples.

        Returns:
            dict[str, float]: A dictionary mapping each sample accession to its normalized
                              median intensity value.
        """
        med_map: dict[str, float] = {}
        for _, batch_df in self.iter_samples(1000, ["sample_accession", "intensity"]):
            meds = batch_df.groupby(["sample_accession"])["intensity"].median()
            med_map.update(meds.to_dict())
        global_med = np.median([med for med in med_map.values()])
        for sample, med in med_map.items():
            med_map[sample] = med / global_med
        return med_map

    def get_report_condition_from_database(self, cons: list, columns: list = None) -> pd.DataFrame:
        """
        Retrieves a standardized report from the database for specified conditions.

        Parameters:
            cons (list): A list of conditions to filter the report.
            columns (list, optional): A list of column names to include in the report.
                                      If None, all columns are included.

        Returns:
            pd.DataFrame: A DataFrame containing the report with standardized column names.
        """
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            f"""SELECT {cols} FROM parquet_db WHERE condition IN {tuple(cons)}"""
        )
        report = database.df()
        return Feature.standardize_df(report)

    def iter_conditions(
        self, conditions: int = 10, columns: list = None
    ) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """
        Iterates over experimental conditions in batches, yielding each batch along with its
        corresponding report DataFrame.

        Parameters:
            conditions (int, optional): The number of conditions to include in each batch. Defaults to 10.
            columns (list, optional): A list of column names to include in the report. If None, all columns are included.

        Yields:
            Iterator[tuple[list[str], pd.DataFrame]]: An iterator over tuples, each containing a list of condition names
                                                      and a DataFrame with the report for those conditions.
        """
        condition_list = self.get_unique_conditions()
        ref_list = [
            condition_list[i : i + conditions] for i in range(0, len(condition_list), conditions)
        ]
        for refs in ref_list:
            batch_df = self.get_report_condition_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_conditions(self) -> list[str]:
        """
        Retrieves a list of unique experimental conditions from the Parquet database.

        Returns:
            list[str]: A list of unique condition identifiers.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT condition FROM parquet_db").df()
        return unique["condition"].tolist()

    def get_median_map_to_condition(self) -> dict[str, dict[str, float]]:
        """
        Computes a median intensity map for each experimental condition, normalizing
        the median intensity of each sample within a condition by the mean median
        intensity across all samples in that condition.

        Returns:
            dict[str, dict[str, float]]: A dictionary mapping each condition to another
                                         dictionary, which maps each sample accession
                                         to its normalized median intensity value.
        """
        med_map = {}
        for cons, batch_df in self.iter_conditions(
            1000, ["condition", "sample_accession", "intensity"]
        ):
            for con in cons:
                meds = (
                    batch_df[batch_df["condition"] == con]
                    .groupby(["sample_accession"])["intensity"]
                    .median()
                )
                meds = meds / meds.mean()
                med_map[con] = meds.to_dict()
        return med_map


@log_execution_time(logger)
def peptide_normalization(
    parquet: str,
    sdrf: str,
    min_aa: int,
    min_unique: int,
    remove_ids: str,
    remove_decoy_contaminants: bool,
    remove_low_frequency_peptides: bool,
    output: str,
    skip_normalization: bool,
    nmethod: str,
    pnmethod: str,
    log2: bool,
    save_parquet: bool,
) -> None:
    """
    Perform peptide normalization on a proteomics dataset.

    This function processes a dataset stored in a Parquet file, applying various
    normalization and filtering steps to prepare peptide data for analysis. It
    supports removing contaminants, low-frequency peptides, and user-specified
    proteins, as well as normalizing intensities across technical repetitions and
    conditions. The results can be saved incrementally to CSV and Parquet files.

    Parameters
    ----------
    parquet : str
        Path to the Parquet file containing the dataset.
    sdrf : str
        Path to the SDRF file for quantification details.
    min_aa : int
        Minimum number of amino acids required for peptides.
    min_unique : int
        Minimum number of unique peptides per protein.
    remove_ids : str
        Path to a file with protein IDs to remove.
    remove_decoy_contaminants : bool
        Whether to remove decoys and contaminants.
    remove_low_frequency_peptides : bool
        Whether to remove low-frequency peptides.
    output : str
        Path to the output file for saving results.
    skip_normalization : bool
        Whether to skip normalization steps.
    nmethod : str
        Method for feature-level normalization.
    pnmethod : str
        Method for peptide-level normalization.
    log2 : bool
        Whether to apply log2 transformation to intensities.
    save_parquet : bool
        Whether to save results in Parquet format.
    """

    if os.path.exists(output):
        raise FileExistsError("The output file already exists.")

    if parquet is None:
        raise FileNotFoundError("The file does not exist.")

    feature_normalization = FeatureNormalizationMethod.from_str(nmethod)
    peptide_normalized = PeptideNormalizationMethod.from_str(pnmethod)

    logger.info("Loading data from %s...", parquet)
    feature = Feature(parquet)

    if sdrf:
        technical_repetitions, label, sample_names, choice = analyse_sdrf(sdrf)
    else:
        technical_repetitions, label, sample_names, choice = feature.experimental_inference

    if remove_low_frequency_peptides:
        low_frequency_peptides = feature.low_frequency_peptides

    med_map = {}
    if not skip_normalization and peptide_normalized == PeptideNormalizationMethod.GlobalMedian:
        med_map = feature.get_median_map()
    elif (
        not skip_normalization and peptide_normalized == PeptideNormalizationMethod.ConditionMedian
    ):
        med_map = feature.get_median_map_to_condition()

    # Incremental CSV writing
    write_csv = True
    if write_csv:
        write_csv_task = WriteCSVTask(output)
        write_csv_task.start()

    # Incremental Parquet writing
    if save_parquet:
        writer_parquet_task = WriteParquetTask(output)
        writer_parquet_task.start()

    for samples, df in feature.iter_samples():
        df.dropna(subset=["pg_accessions"], inplace=True)
        for sample in samples:
            # Perform data preprocessing on every sample
            logger.info("%s: Data preprocessing...", str(sample).upper())
            dataset_df = df[df["sample_accession"] == sample].copy()

            # Step1: Parse the identifier of proteins and retain only unique peptides.
            dataset_df = dataset_df[dataset_df["unique"] == 1]
            dataset_df = dataset_df[PARQUET_COLUMNS]

            dataset_df = reformat_quantms_feature_table_quant_labels(dataset_df, label, choice)

            # Step2: Remove lines where intensity or study condition is empty.
            # Step3: Filter peptides with less amino acids than min_aa.
            dataset_df = apply_initial_filtering(dataset_df, min_aa)

            # Step4: Delete low-confidence proteins.
            dataset_df = dataset_df.groupby(PROTEIN_NAME).filter(
                lambda x: len(set(x[PEPTIDE_CANONICAL])) >= min_unique
            )

            # Step5: Filter decoy, contaminants, entrapment
            if remove_decoy_contaminants:
                dataset_df = remove_contaminants_entrapments_decoys(dataset_df)

            # Step6: Filter user-specified proteins
            if remove_ids is not None:
                dataset_df = remove_protein_by_ids(dataset_df, remove_ids)
            dataset_df.rename(columns={INTENSITY: NORM_INTENSITY}, inplace=True)

            # Step7: Normalize at feature level between ms runs (technical repetitions)
            if (
                not skip_normalization
                and nmethod not in ("none", None)
                and technical_repetitions > 1
            ):
                start_time = time.time()
                logger.info(
                    "%s: Normalizing intensities of features using method %s...",
                    str(sample).upper(),
                    nmethod,
                )
                dataset_df = feature_normalization(dataset_df, technical_repetitions)
                # dataset_df = normalize_runs(dataset_df, technical_repetitions, nmethod)
                elapsed = time.time() - start_time
                logger.info(
                    "%s: Number of features after normalization: %d (completed in %.2f seconds)",
                    str(sample).upper(),
                    len(dataset_df.index),
                    elapsed,
                )
            # Step8: Merge peptidoforms across fractions and technical repetitions
            dataset_df = get_peptidoform_normalize_intensities(dataset_df)
            logger.info(
                "%s: Number of peptides after peptidoform selection: %d",
                str(sample).upper(),
                len(dataset_df.index),
            )

            if len(dataset_df[FRACTION].unique().tolist()) > 1:
                start_time = time.time()
                logger.info("%s: Merging features across fractions...", str(sample).upper())
                dataset_df = merge_fractions(dataset_df)
                elapsed = time.time() - start_time
                logger.info(
                    "%s: Number of features after merging fractions: %d (completed in %.2f seconds)",
                    str(sample).upper(),
                    len(dataset_df.index),
                    elapsed,
                )
            # Step9: Normalize the data.
            if not skip_normalization:
                dataset_df = peptide_normalized(dataset_df, sample, med_map)

            # Step10: Remove peptides with low frequency.
            if remove_low_frequency_peptides and len(sample_names) > 1:
                dataset_df.set_index([PROTEIN_NAME, PEPTIDE_CANONICAL], drop=True, inplace=True)
                dataset_df = dataset_df[
                    ~dataset_df.index.isin(low_frequency_peptides)
                ].reset_index()
                logger.info(
                    "%s: Peptides after removing low frequency peptides: %d",
                    str(sample).upper(),
                    len(dataset_df.index),
                )

            # Step11: Assembly peptidoforms to peptides.
            start_time = time.time()
            logger.info("%s: Summing all peptidoforms per sample...", str(sample).upper())
            dataset_df = sum_peptidoform_intensities(dataset_df)
            elapsed = time.time() - start_time
            logger.info(
                "%s: Number of peptides after selection: %d (completed in %.2f seconds)",
                str(sample).upper(),
                len(dataset_df.index),
                elapsed,
            )
            # Step12: Intensity transformation to log.
            if log2:
                dataset_df[NORM_INTENSITY] = np.log2(dataset_df[NORM_INTENSITY])

            logger.info("%s: Saving the normalized peptide intensities...", str(sample).upper())

            if save_parquet:
                writer_parquet_task.write(dataset_df)
            if write_csv:
                write_csv_task.write(dataset_df)

    if write_csv:
        write_csv_task.close()
    if save_parquet:
        writer_parquet_task.close()
