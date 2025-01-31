import os
import re

from typing import Iterator, Optional

import pandas as pd
import numpy as np
import duckdb

from ibaqpy.model.quantification_type import (QuantificationCategory, IsobaricLabel)
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


def parse_uniprot_accession(uniprot_id: str) -> str:
    """
    Parse the uniprot accession from the uniprot id in the form of
    tr|CONTAMINANT_Q3SX28|CONTAMINANT_TPM2_BOVIN and convert to CONTAMINANT_Q3SX28
    :param uniprot_id: uniprot id
    :return: uniprot accession
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
    This function returns a peptide sequence without the modification information
    :param peptide_sequence: peptide sequence with mods
    :return: peptide sequence
    """
    clean_peptide = re.sub(r"[\(\[].*?[\)\]]", "", peptide_sequence)
    clean_peptide = clean_peptide.replace(".", "").replace("-", "")
    return clean_peptide


def analyse_sdrf(sdrf_path: str) -> tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
    """
    This function is aimed to parse SDRF and return four objects:
    1. sdrf_df: A dataframe with channels and references annoted.
    2. label: Label type of the experiment. LFQ, TMT or iTRAQ.
    3. sample_names: A list contains all sample names.
    4. choice: A dictionary caontains key-values between channel
        names and numbers.
    :param sdrf_path: File path of SDRF.
    :return:
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
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants = ["CONTAMINANT", "ENTRAP", "DECOY"]
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex)]


def remove_protein_by_ids(
    dataset: pd.DataFrame, protein_file: str, protein_field=PROTEIN_NAME
) -> pd.DataFrame:
    """
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param protein_file: contaminants file
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants_reader = open(protein_file, "r")
    contaminants = contaminants_reader.read().split("\n")
    contaminants = [cont for cont in contaminants if cont.strip()]
    cregex = "|".join(contaminants)
    return dataset[~dataset[protein_field].str.contains(cregex, regex=True)]


def reformat_quantms_feature_table_quant_labels(data_df: pd.DataFrame, label: QuantificationCategory, choice: Optional[IsobaricLabel]) -> pd.DataFrame:
    """
    Reformat (a subset of) a ``quantms`` feature table for consistent processing.

    :param data_df: Feature data in dataframe.
    :param label: Label type of the experiment.
    :param choice: Choice dict for a label type.
    :return: Processed data.
    """

    data_df = data_df.rename(columns=parquet_map)
    data_df[PROTEIN_NAME] = data_df[PROTEIN_NAME].str.join(";")
    if label == QuantificationCategory.LFQ:
        data_df.drop(CHANNEL, inplace=True, axis=1)
    else:
        data_df[CHANNEL] = data_df[CHANNEL].map(choice.channels())

    return data_df


def apply_initial_filtering(data_df: pd.DataFrame, min_aa: int) -> pd.DataFrame:
    # Remove 0 intensity signals from the data
    data_df = data_df[data_df[INTENSITY] > 0]

    data_df = data_df[
        (data_df["Condition"] != "Empty") | (data_df["Condition"].isnull())
    ]

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
    Merge features across fractions.
    :param dataset: dataset including all properties
    :return:
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
    Select the best peptidoform for the same sample and the same replicates. A peptidoform is the combination of
    a (PeptideSequence + Modifications) + Charge state.
    :param dataset: dataset including all properties
    :param higher_intensity: select based on normalize intensity, if false based on best scored peptide
    :return:
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
    Sum the peptidoform intensities for all peptidofrom across replicates of the same sample.
    :param dataset: Dataframe to be analyzed
    :return: dataframe with the intensities
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

    def standardize_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.rename({"protein_accessions": "pg_accessions", "charge": "precursor_charge"}, axis=1)

    @property
    def experimental_inference(self) -> tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
        self.labels = self.get_unique_labels()
        self.label, self.choice = QuantificationCategory.classify(self.labels)
        self.technical_repetitions = self.get_unique_tec_reps()
        return len(self.technical_repetitions), self.label, self.samples, self.choice

    @property
    def low_frequency_peptides(self, percentage=0.2) -> tuple:
        """Return peptides with low frequency"""
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
            print(e)
            exit("Some errors occurred when parsing pg_accessions column in feature parquet!")
        f_table.set_index(["sequence", "pg_accessions"], inplace=True)
        f_table.drop(
            f_table[f_table["count"] >= (percentage * len(self.samples))].index,
            inplace=True,
        )
        f_table.reset_index(inplace=True)
        return tuple(zip(f_table["pg_accessions"], f_table["sequence"]))

    @staticmethod
    def csv2parquet(csv):
        parquet_path = os.path.splitext(csv)[0] + ".parquet"
        duckdb.read_csv(csv).to_parquet(parquet_path)

    def get_report_from_database(self, samples: list, columns: list = None):
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param columns: A list of columns
        :param samples: A list of samples
        :return: The report
        """
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            """SELECT {} FROM parquet_db WHERE sample_accession IN {}""".format(
                cols, tuple(samples)
            )
        )
        report = database.df()
        return self.standardize_df(report)

    def iter_samples(self, sample_num: int = 20, columns: list = None) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """
        :params sample_num: The number of samples being processed at the same time (default 20)
        :yield: _description_
        """
        ref_list = [self.samples[i : i + sample_num] for i in range(0, len(self.samples), sample_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_samples(self) -> list[str]:
        """
        return: A list of samples.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT sample_accession FROM parquet_db").df()
        return unique["sample_accession"].tolist()

    def get_unique_labels(self) -> list[str]:
        """
        return: A list of labels.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT channel FROM parquet_db").df()
        return unique["channel"].tolist()

    def get_unique_tec_reps(self) -> list[int]:
        """
        return: A list of labels.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT run FROM parquet_db").df()
        try:
            if unique["run"].str.contains('_').all():
                unique["run"] = unique["run"].str.split("_").str.get(1)
                unique["run"] = unique["run"].astype("int")
            else:
                unique["run"] = unique["run"].astype("int")
        except ValueError as e:
            raise ValueError(f"Some errors occurred when getting technical repetitions: {e}") from e

        return unique["run"].tolist()

    def get_median_map(self) -> dict[str, float]:
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
        This function loads the report from the duckdb database for a group of ms_runs.
        :param columns: A list of columns
        :param cons: A list of conditions in
        :return: The report
        """
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            f"""SELECT {cols} FROM parquet_db WHERE condition IN {tuple(cons)}"""
        )
        report = database.df()
        return self.standardize_df(report)

    def iter_conditions(self, conditions: int = 10, columns: list = None) -> Iterator[tuple[list[str], pd.DataFrame]]:
        condition_list = self.get_unique_conditions()
        ref_list = [
            condition_list[i : i + conditions] for i in range(0, len(condition_list), conditions)
        ]
        for refs in ref_list:
            batch_df = self.get_report_condition_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_conditions(self) -> list[str]:
        """
        return: A list of conditions.
        """
        unique = self.parquet_db.sql("SELECT DISTINCT condition FROM parquet_db").df()
        return unique["condition"].tolist()

    def get_median_map_to_condition(self) -> dict[str, dict[str, float]]:
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

    :param parquet: Parquet file with features
    :param sdrf: SDRF file
    :param min_aa: Min amino acids
    :param min_unique: Min of unique peptides
    :param remove_ids: Remove features for the given proteins
    :param remove_decoy_contaminants: Remove contaminants and entrapment
    :param remove_low_frequency_peptides: Remove low frecuency peptides
    :param output: Output file
    :param skip_normalization: Skip normalization
    :param nmethod: normalization method for features
    :param pnmethod: peptide normalization method
    :param log2: log intensities for features before normalizing
    :param save_parquet: Save to parque file.
    """
    if os.path.exists(output):
        raise FileExistsError("The output file already exists.")

    if parquet is None:
        raise FileNotFoundError("The file does not exist.")

    feature_normalization = FeatureNormalizationMethod.from_str(nmethod)
    peptide_normalization = PeptideNormalizationMethod.from_str(pnmethod)

    print("Loading data..")
    feature = Feature(parquet)

    if sdrf:
        technical_repetitions, label, sample_names, choice = analyse_sdrf(sdrf)
    else:
        technical_repetitions, label, sample_names, choice = feature.experimental_inference

    if remove_low_frequency_peptides:
        low_frequency_peptides = feature.low_frequency_peptides

    med_map = {}
    if not skip_normalization and peptide_normalization == PeptideNormalizationMethod.GlobalMedian:
        med_map = feature.get_median_map()
    elif not skip_normalization and peptide_normalization == PeptideNormalizationMethod.ConditionMedian:
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
            print(f"{str(sample).upper()}: Data preprocessing...")
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
            if not skip_normalization and nmethod not in ("none", None) and technical_repetitions > 1:
                print(f"{str(sample).upper()}: Normalize intensities of features.. ")
                dataset_df = feature_normalization(dataset_df, technical_repetitions)
                # dataset_df = normalize_runs(dataset_df, technical_repetitions, nmethod)
                print(
                    f"{str(sample).upper()}: Number of features after normalization: {len(dataset_df.index)}"
                )
            # Step8: Merge peptidoforms across fractions and technical repetitions
            dataset_df = get_peptidoform_normalize_intensities(dataset_df)
            print(
                f"{str(sample).upper()}: Number of peptides after peptidofrom selection: {len(dataset_df.index)}"
            )

            if len(dataset_df[FRACTION].unique().tolist()) > 1:
                print(f"{str(sample).upper()}: Merge features across fractions.. ")
                dataset_df = merge_fractions(dataset_df)
                print(
                    f"{str(sample).upper()}: Number of features after merging fractions: {len(dataset_df.index)}"
                )
            # Step9: Normalize the data.
            if not skip_normalization:
                dataset_df = peptide_normalization(
                    dataset_df,
                    sample,
                    med_map
                )


            # Step10: Remove peptides with low frequency.
            if remove_low_frequency_peptides and len(sample_names) > 1:
                dataset_df.set_index([PROTEIN_NAME, PEPTIDE_CANONICAL], drop=True, inplace=True)
                dataset_df = dataset_df[
                    ~dataset_df.index.isin(low_frequency_peptides)
                ].reset_index()
                print(
                    f"{str(sample).upper()}: Peptides after remove low frequency peptides: {len(dataset_df.index)}"
                )

            # Step11: Assembly peptidoforms to peptides.
            print(f"{str(sample).upper()}: Sum all peptidoforms per Sample...")
            dataset_df = sum_peptidoform_intensities(dataset_df)
            print(
                f"{str(sample).upper()}: Number of peptides after selection: {len(dataset_df.index)}"
            )
            # Step12: Intensity transformation to log.
            if log2:
                dataset_df[NORM_INTENSITY] = np.log2(dataset_df[NORM_INTENSITY])

            print(f"{str(sample).upper()}: Save the normalized peptide intensities...")

            if save_parquet:
                writer_parquet_task.write(dataset_df)
            if write_csv:
                write_csv_task.write(dataset_df)

    if write_csv:
        write_csv_task.close()
    if save_parquet:
        writer_parquet_task.close()
