#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import click
import numpy as np
import pandas as pd
import os
from pandas import DataFrame

from ibaq.ibaqpy_commons import *


def parse_uniprot_accession(uniprot_id: str) -> str:
    """
    Parse the uniprot accession from the uniprot id in the form of
    tr|CONTAMINANT_Q3SX28|CONTAMINANT_TPM2_BOVIN and convert to CONTAMINANT_TPM2_BOVIN
    :param uniprot_id: uniprot id
    :return: uniprot accession
    """
    uniprot_list = uniprot_id.split(";")
    result_uniprot_list = []
    for accession in uniprot_list:
        if accession.count("|") == 2:
            accession = accession.split("|")[2]
        result_uniprot_list.append(accession)
    return ";".join(result_uniprot_list)


def print_dataset_size(dataset: DataFrame, message: str, verbose: bool) -> None:
    if verbose:
        print(message + str(len(dataset.index)))


def print_help_msg(command) -> None:
    """
    Print help information
    :param command: command to print helps
    :return: print help
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))



def get_peptidoform_normalize_intensities(dataset: DataFrame, higher_intensity: bool = True) -> DataFrame:
    """
    Select the best peptidoform for the same sample and the same replicates. A peptidoform is the combination of
    a (PeptideSequence + Modifications) + Charge state.
    :param dataset: dataset including all properties
    :param higher_intensity: select based on normalize intensity, if false based on best scored peptide
    :return:
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    if higher_intensity:
        dataset = dataset.loc[dataset.groupby([PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE])[
            NORM_INTENSITY].idxmax()].reset_index(drop=True)
    else:
        dataset = dataset.loc[dataset.groupby([PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, CONDITION, BIOREPLICATE])[
            SEARCH_ENGINE].idxmax()].reset_index(drop=True)
    print(dataset)
    return dataset


def sum_peptidoform_intensities(dataset: DataFrame) -> DataFrame:
    """
    Sum the peptidoform intensities for all peptidofrom across replicates of the same sample.
    :param dataset: Dataframe to be analyzed
    :return: dataframe with the intensities
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    normalize_df = dataset.groupby([PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION])[NORM_INTENSITY].sum()
    normalize_df = normalize_df.reset_index()
    normalize_df = pd.merge(normalize_df,
                            dataset[[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION]], how='left',
                            on=[PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE, CONDITION])

    return normalize_df


def average_peptide_intensities(dataset: DataFrame) -> DataFrame:
    """
    Median the intensities of all the peptidoforms for a specific peptide sample combination.
    :param dataset: Dataframe containing all the peptidoforms
    :return: New dataframe
    """
    dataset_df = dataset.groupby([PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION])[NORM_INTENSITY].median()
    dataset_df = dataset_df.reset_index()
    dataset_df = pd.merge(dataset_df, dataset[[PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION]], how='left',
                          on=[PEPTIDE_CANONICAL, SAMPLE_ID, CONDITION])

    return dataset_df


def peptide_intensity_normalization(dataset_df: DataFrame, field: str, class_field: str, scaling_method: str):
    """
    Normalize the peptide intensities using different methods.
    :param dataset_df: dataframe with the data
    :param field: field to normalize
    :param class_field: field to use as class
    :param scaling_method: method to use for the normalization
    :return:
    """
    if scaling_method == 'qnorm':
        # pivot to have one col per sample
        normalize_df = pd.pivot_table(dataset_df, index=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION],
                                      columns=class_field, values=field, aggfunc={field: np.mean})
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION])
        normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
        normalize_df = normalize_df[normalize_df[NORM_INTENSITY].notna()]
        return normalize_df

    return dataset_df


def impute_peptide_intensities(dataset_df, field, class_field):
    """
    Impute the missing values using different methods.
    :param dataset_df: dataframe with the data
    :param field: field to impute
    :param class_field: field to use as class
    :return:
    """
    normalize_df = pd.DataFrame()
    # group by condition to detect missing values
    for c, g in dataset_df.groupby(CONDITION):
        # pivot to have one col per sample
        group_normalize_df = pd.pivot_table(g, index=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION],
                                            columns=class_field, values=field, aggfunc={field: np.mean})

        # no missing values group -> only one sample
        if len(group_normalize_df.columns) < 2:
            group_normalize_df = group_normalize_df.reset_index()
            group_normalize_df = group_normalize_df.melt(id_vars=[PEPTIDE_CANONICAL, PROTEIN_NAME, CONDITION])
            group_normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
            normalize_df = normalize_df.append(group_normalize_df, ignore_index=True)

    return normalize_df


def remove_extension_file(filename: str) -> str:
    """
  The filename can have
  :param filename:
  :return:
  """
    return filename.replace('.raw', '').replace('.RAW', '').replace('.mzML', '').replace('.wiff', '')


def get_study_accession(sample_id: str) -> str:
    """
  Get the project accession from the Sample accession. The function expected a sample accession in the following
  format PROJECT-SAMPLEID
  :param sample_id: Sample Accession
  :return: study accession
  """
    return sample_id.split('-')[0]


def get_reference_name(reference_spectrum: str) -> str:
    """
    Get the reference name from Reference column. The function expected a reference name in the following format eg.
    20150820_Haura-Pilot-TMT1-bRPLC03-2.mzML_controllerType=0 controllerNumber=1 scan=16340
    :param reference_spectrum:
    :return: reference name
    """
    return re.split(r'\.mzML|\.MZML|\.raw|\.RAW', reference_spectrum)[0]


@click.command()
@click.option("-m", "--msstats", help="MsStats file import generated by quantms", required=True)
@click.option("-s", "--sdrf", help="SDRF file import generated by quantms", required=True)
@click.option("--compress", help="Read all files compress", is_flag=True)
@click.option("--chunksize", help="The number of rows of MSstats read using pandas streaming", default=10*1024*1024)
@click.option("--min_aa", help="Minimum number of amino acids to filter peptides", default=7)
@click.option("--min_unique", help="Minimum number of unique peptides to filter proteins", default=2)
@click.option("--contaminants", help="Contaminants and high abundant proteins to be removed")
@click.option("--skip_normalization", help="Skip normalization step", is_flag=True, default=False)
@click.option('--nmethod', help="Normalization method used to normalize intensities for all samples (options: qnorm)",
              default="qnorm")
@click.option("--impute", help="Impute the missing values using MissForest", is_flag=True)
@click.option("--pnormalization", help="Normalize the peptide intensities using different methods (options: qnorm)",
              is_flag=True)
@click.option("--output", help="Peptide intensity file including other all properties for normalization")
@click.option("--compress", help="Read the input peptides file in compress gzip file", is_flag=True)
@click.option("--log2", help="Transform to log2 the peptide intensity values before normalization", is_flag=True)

def peptide_normalization(msstats: str, sdrf: str, contaminants: str, output: str, skip_normalization: bool,
                          min_aa: int, min_unique: int, nmethod: str, impute: bool, pnormalization: bool,
                          compress: bool, chunksize: int, log2: bool) -> None:
    """Intensity normalization of stream processing peptide performed from MSstats

    :param msstats:
    :param sdrf:
    :param contaminants:
    :param output:
    :param skip_normalization:
    :param min_aa:
    :param min_unique:
    :param nmethod:
    :param impute:
    :param pnormalization:
    :param compress:
    :param chunksize:
    :param log2:
    :return:
    """

    if msstats is None or output is None:
        print_help_msg(peptide_normalization)
        exit(1)

    pd.set_option('display.max_columns', None)
    print("Loading data..")
    compression_method = 'gzip' if compress else None

    # Read the sdrf file
    sdrf_df = pd.read_csv(sdrf, sep='\t', compression=compression_method)
    sdrf_df[REFERENCE] = sdrf_df['comment[data file]'].apply(remove_extension_file)
    print(sdrf_df)

    # Determine label type
    labels = set(sdrf_df['comment[label]'])
    with open(msstats, 'r') as file:
        header_row = file.readline().strip().split(',')  
        msstats_columns = header_row[1:]
    if CHANNEL not in msstats_columns:
        label = "LFQ"
    elif 'TMT' in ','.join(labels) or 'tmt' in ','.join(labels):
        if (len(labels) > 11 or "TMT134N" in labels or "TMT133C" in labels
                or "TMT133N" in labels or "TMT132C" in labels or "TMT132N" in labels):
            choice = TMT16plex
        elif len(labels) == 11 or "TMT131C" in labels:
            choice = TMT11plex
        elif len(labels) > 6:
            choice = TMT10plex
        else:
            choice = TMT6plex
        choice = pd.DataFrame.from_dict(choice, orient='index', columns=[CHANNEL]).reset_index().rename(
            columns={'index': 'comment[label]'})
        sdrf_df = sdrf_df.merge(choice, on='comment[label]', how='left')
        label = "TMT"
    elif 'ITRAQ' in ','.join(labels) or 'itraq' in ','.join(labels):
        if len(labels) > 4:
            choice = ITRAQ8plex
        else:
            choice = ITRAQ4plex
        choice = pd.DataFrame.from_dict(choice, orient='index', columns=[CHANNEL]).reset_index().rename(
            columns={'index': 'comment[label]'})
        sdrf_df = sdrf_df.merge(choice, on='comment[label]', how='left')
        label = "ITRAQ"
    else:
        print("Warning: Only support label free, TMT and ITRAQ experiment!")
        exit(1)
    sample_names = sdrf_df["source name"].unique().tolist()

    ## TODO: Stream processing to obtain strong proteins with more than 2 uniqe peptides
    print("IBAQPY WARNING: Writing files into ibaqpyTemp...")
    if not os.path.exists("ibaqpyTemp/"):
        os.mkdir("ibaqpyTemp/")
    msstats_chunks = pd.read_csv(msstats, sep=',', compression=compression_method, dtype={CONDITION: 'category',
                                 ISOTOPE_LABEL_TYPE: 'category'}, chunksize=chunksize)
    unique_peptides = {}
    canonical_dict = {}
    group_intensities = {}
    quantile = {}
    for msstats_df in msstats_chunks:
        msstats_df.rename(
        columns={'ProteinName': PROTEIN_NAME, 'PeptideSequence': PEPTIDE_SEQUENCE, 'PrecursorCharge': PEPTIDE_CHARGE,
                 'Run': RUN,
                 'Condition': CONDITION, 'Intensity': INTENSITY}, inplace=True)
        msstats_df = msstats_df[msstats_df[INTENSITY] > 0]
        if PEPTIDE_CANONICAL not in msstats_df.columns:
            modified_seqs = msstats_df[PEPTIDE_SEQUENCE].unique().tolist()
            canonical_seqs = [get_canonical_peptide(i) for i in modified_seqs]
            inner_canonical_dict = dict(zip(modified_seqs, canonical_seqs))
            canonical_dict.update(inner_canonical_dict)
            msstats_df[PEPTIDE_CANONICAL] = msstats_df.apply(lambda x: inner_canonical_dict[x[PEPTIDE_SEQUENCE]], axis=1)
        # Filter peptides with less amino acids than min_aa (default: 7)
        msstats_df = msstats_df[msstats_df.apply(lambda x: len(x[PEPTIDE_CANONICAL]) >= min_aa, axis=1)]
        msstats_df[PROTEIN_NAME] = msstats_df[PROTEIN_NAME].apply(parse_uniprot_accession)

        if FRACTION not in msstats_df.columns:
            msstats_df[FRACTION] = 1
            msstats_df = msstats_df[
                [PROTEIN_NAME, PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, INTENSITY, REFERENCE, CONDITION, RUN,
                BIOREPLICATE, FRACTION, FRAGMENT_ION, ISOTOPE_LABEL_TYPE]]

        # Merged the SDRF with Resulted file
        if label == "LFQ":
            msstats_df[REFERENCE] = msstats_df[REFERENCE].apply(remove_extension_file)
            result_df = pd.merge(msstats_df, sdrf_df[['source name', REFERENCE]], how='left', on=[REFERENCE])
        elif label == "TMT":
            msstats_df[REFERENCE] = msstats_df[REFERENCE].apply(get_reference_name)
            result_df = pd.merge(msstats_df, sdrf_df[['source name', REFERENCE, CHANNEL]], how='left',
                                on=[REFERENCE, CHANNEL])
            result_df = result_df[result_df["Condition"] != "Empty"]
            result_df.rename(columns={'Charge': PEPTIDE_CHARGE}, inplace=True)
        elif label == 'ITRAQ':
            msstats_df[REFERENCE] = msstats_df[REFERENCE].apply(get_reference_name)
            result_df = pd.merge(msstats_df, sdrf_df[['source name', REFERENCE, CHANNEL]], how='left',
                                on=[REFERENCE, CHANNEL])
            result_df = result_df[result_df["Condition"] != "Empty"]
            result_df.rename(columns={'Charge': PEPTIDE_CHARGE}, inplace=True)

        result_df.rename(columns={'source name': SAMPLE_ID}, inplace=True)
        result_df[STUDY_ID] = result_df[SAMPLE_ID].apply(get_study_accession)

        # Write CSVs by Sample ID
        for sample in sample_names:
            file_name = f"ibaqpyTemp/{sample}.csv"
            write_mode = "a" if os.path.exists(file_name) else "w"
            header = False if os.path.exists(file_name) else True
            result_df[result_df[SAMPLE_ID] == sample].to_csv(file_name, index=False, header = header, mode=write_mode)
        unique_df = result_df.groupby(PEPTIDE_CANONICAL).filter(lambda x: len(set(x[PROTEIN_NAME])) == 1)[[
                              PEPTIDE_CANONICAL, PROTEIN_NAME]]
        unique_dict = dict(zip(unique_df[PEPTIDE_CANONICAL], unique_df[PROTEIN_NAME]))
        for i in unique_dict.keys():
            if i in unique_peptides.keys() and unique_dict[i] != unique_peptides[i]:
                unique_peptides.pop(i)
            else:
                unique_peptides[i] = unique_dict[i]

    proteins_list = list(unique_peptides.values())
    count_dict = {element: proteins_list.count(element) for element in set(proteins_list)}
    strong_proteins = [element for element in count_dict if count_dict[element] >= min_unique]
    del proteins_list, count_dict
    print(f"Number of unique peptides: {len(list(unique_peptides.keys()))}")
    print(f"Number of strong proteins: {len(strong_proteins)}")

    ## TODO: Filter proteins with less unique peptides than min_unique (default: 2)
    for sample in sample_names:
        print(f"{sample} -> Filter out proteins containing unique peptides fewer than {min_unique}..")
        msstats_df = pd.read_csv(f"ibaqpyTemp/{sample}.csv", sep=',')
        msstats_df = msstats_df[msstats_df[PROTEIN_NAME].isin(strong_proteins)]
        print(f"{sample} -> Logarithmic if specified..")
        msstats_df.loc[msstats_df.Intensity == 0, INTENSITY] = 1
        msstats_df[NORM_INTENSITY] = np.log2(msstats_df[INTENSITY]) if log2 else msstats_df[INTENSITY]
        msstats_df.to_csv(f"ibaqpyTemp/{sample}.csv", index=False, sep=',')
        if not skip_normalization:
            if nmethod == 'msstats':
                if label in ["TMT", "ITRAQ"]:
                    g = msstats_df.groupby(['Run', 'Channel'])
                else:
                    g = msstats_df.groupby(['Run', 'Fraction'])
                for name, group in g:
                    group_intensity = group[NORM_INTENSITY].tolist()
                    if name not in group_intensities:
                        group_intensities[name] = group_intensity
                    else:
                        group_intensities.update({name: group_intensities[NORM_INTENSITY] + group_intensity})
            elif nmethod == 'qnorm':

                def recalculate(original, count, value):
                    return (original*count+value)/(count+1), count + 1
                
                dic = msstats_df[NORM_INTENSITY].dropna().sort_values(ascending=False).reset_index(drop=True).to_dict()
                if len(quantile) == 0:
                    quantile = {k: (v, 1) for k, v in dic.items()}
                else:
                    update = min(len(quantile), len(dic))
                    for i in range(0, update):
                        original, count = quantile[i]
                        quantile[i] = recalculate(original, count, dic[i])
                    if len(dic) <= len(quantile):
                        continue
                    else:
                        quantile.update({k: (v, 1) for k, v in dic.items() if k >= update})

    def normalization(dataset_df, label, sample, skip_normalization, nmethod, pnormalization):
        # Remove high abundant and contaminants proteins and the outliers
        if contaminants is not None:
            print(f"{sample} -> Remove contaminants...")
            dataset_df = remove_contaminants_decoys(dataset_df, contaminants)
            print(f"{sample} -> Peptides after contaminants removal: {len(dataset_df[PEPTIDE_SEQUENCE].unique().tolist())}")

        if not skip_normalization:
            print(f"{sample} -> Normalize intensities.. ")
            class_field = SAMPLE_ID
            field = NORM_INTENSITY
            if nmethod == 'msstats':
                # For ISO normalization
                if label in ["TMT", "ITRAQ"]:
                    dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(lambda x: x[field] - group_intensities[(x["Run"], x["Channel"])] +
                                                                         median_baseline, axis = 1)
                else:
                    dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(lambda x: x[field] - group_intensities[(x["Run"], x["Fraction"])] +
                                                                         np.median([group_intensities[i] for i in group_intensities.keys()
                                                                         if i[1] == x["Fraction"]]), axis = 1)
            elif nmethod == 'qnorm':
                # pivot to have one col per sample
                ref_dict = dataset_df[NORM_INTENSITY].dropna().drop_duplicates().sort_values(ascending=False).reset_index(drop=True).to_dict()
                ref_dict = {v: norm_intensity[k] for k, v in ref_dict.items()}
                dataset_df.loc[:, NORM_INTENSITY] = dataset_df.apply(lambda x: ref_dict[x[NORM_INTENSITY]] if x[NORM_INTENSITY] in ref_dict.keys()
                                                                     else np.nan, axis = 1)

        print(f"{sample} -> Select the best peptidoform across fractions...")
        print(f"{sample} -> Number of peptides before peptidofrom selection: " + str(len(dataset_df.index)))
        dataset_df = get_peptidoform_normalize_intensities(dataset_df)
        print(f"{sample} -> Number of peptides after peptidofrom selection: " + str(len(dataset_df.index)))

        # Add the peptide sequence canonical without the modifications
        if PEPTIDE_CANONICAL not in dataset_df.columns:
            print(f"{sample} -> Add Canonical peptides to the dataframe...")
            dataset_df[PEPTIDE_CANONICAL] = dataset_df[PEPTIDE_SEQUENCE].apply(lambda x: get_canonical_peptide(x))

        print(f"{sample} -> Sum all peptidoforms per Sample...")
        print(f"{sample} -> Number of peptides before sum selection: " + str(len(dataset_df.index)))
        dataset_df = sum_peptidoform_intensities(dataset_df)
        print(f"{sample} -> Number of peptides after sum: " + str(len(dataset_df.index)))

        print(f"{sample} -> Average all peptidoforms per Peptide/Sample...")
        print(f"{sample} -> Number of peptides before average: " + str(len(dataset_df.index)))
        dataset_df = average_peptide_intensities(dataset_df)
        print(f"{sample} -> Number of peptides after average: " + str(len(dataset_df.index)))

        # Perform imputation using Random Forest in Peptide Intensities
        # TODO: Check if this is necessary (Probably we can do some research if imputation at peptide level is necessary
        # if impute:
        #     dataset_df = impute_peptide_intensities(dataset_df, field=NORM_INTENSITY, class_field=SAMPLE_ID)

        if pnormalization:
            print(f"{sample} -> Normalize at Peptide level...")
            dataset_df = peptide_intensity_normalization(dataset_df, field=NORM_INTENSITY, class_field=SAMPLE_ID,
                                                        scaling_method=nmethod)
        dataset_df = dataset_df.drop_duplicates()
        dataset_df = dataset_df[dataset_df[NORM_INTENSITY].notna()]
        return dataset_df

    ## TODO: Peptide intensity normalization
    peptides_count = {}
    if not skip_normalization:
        if nmethod == 'qnorm':
            norm_intensity = {k: v[0] for k, v in quantile.items()}
        elif nmethod == 'msstats':
            # For ISO normalization
            print(f"Label -> {label}")
            if label in ["TMT", "ITRAQ"]:
                median_baseline = np.median(list(set(sum(group_intensities.values(), []))))
                group_intensities = {key: np.median(list(values)) for key, values in group_intensities.items()}
            else:
                fractions = [i[1] for i in group_intensities.keys()]
                fraction_median = {}
                for fraction in fractions:
                    fraction_keys = [i for i in group_intensities.keys() if i[1] == fraction]
                    fraction_intensities = []
                    for key in fraction_keys:
                        fraction_intensities.extend(group_intensities[key])
                    fraction_median[fraction] = np.median(fraction_intensities)
                group_intensities = {key: np.median(values) for key, values in group_intensities.items()}
    for sample in sample_names:
        dataset_df = pd.read_csv(f"ibaqpyTemp/{sample}.csv", sep=',')
        norm_df = normalization(dataset_df, label, sample, skip_normalization, nmethod, pnormalization)
        sample_peptides = norm_df[PEPTIDE_CANONICAL].unique().tolist()
        peptides_count = {peptide: peptides_count.get(peptide, 0) + 1 for peptide in sample_peptides}
        norm_df.to_csv(f"ibaqpyTemp/{sample}.csv", sep = ",", index=False)
    del group_intensities, quantile

    sample_number = len(sample_names)
    peptides_count = {k: v / sample_number for k, v in peptides_count.items() if v / sample_number >= 0.2 and v > 1}
    strong_peptides = list(peptides_count.keys())
    del peptides_count, strong_proteins

    # Filter low-frequency peptides
    print("IBAQPY WARNING: Writing normalized intensities into CSV...")
    for sample in sample_names:
        dataset_df = pd.read_csv(f"ibaqpyTemp/{sample}.csv", sep=',')
        # Filter low-frequency peptides, which indicate whether the peptide occurs less than 20% in all samples or
        # only in one sample
        dataset_df = dataset_df[dataset_df[PEPTIDE_CANONICAL].isin(strong_peptides)]
        dataset_df = dataset_df[[PEPTIDE_CANONICAL, PROTEIN_NAME, SAMPLE_ID, NORM_INTENSITY, CONDITION]]
        write_mode = "a" if os.path.exists(output) else "w"
        header = False if os.path.exists(output) else True
        dataset_df.to_csv(output, index=False, header=header, mode=write_mode)
        dataset_df.to_csv(f"ibaqpyTemp/{sample}.csv", sep=',', index=False)
    del strong_peptides


if __name__ == '__main__':
    peptide_normalization()
