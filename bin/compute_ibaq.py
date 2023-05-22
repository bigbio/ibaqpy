#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math

import click
import pandas as pd
from pandas import DataFrame, Series
from pyopenms import *

from ibaq.ibaqpy_commons import PROTEIN_NAME, IBAQ, IBAQ_LOG, IBAQ_PPB, NORM_INTENSITY, SAMPLE_ID, IBAQ_NORMALIZED, \
    CONDITION
from ibaq.ibaqpy_commons import plot_distributions, plot_box_plot

def print_help_msg(command):
    """
    Print the help of the command
    :param command: command
    :return:
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


def normalize(group):
    group[IBAQ_NORMALIZED] = group[IBAQ] / group[IBAQ].sum()
    return group


def normalize_ibaq(res: DataFrame) -> DataFrame:
    """
    Normalize the ibaq values using the total ibaq of the sample. The resulted
    ibaq values are then multiplied by 100'000'000 (PRIDE database noramalization)
    for the ibaq ppb and log10 shifted by 10 (ProteomicsDB)
    :param res: Dataframe
    :return:
    """

    res = res.groupby([SAMPLE_ID, CONDITION]).apply(normalize)

    # Normalization method used by Proteomics DB 10 + log10(ibaq/sum(ibaq))
    res[IBAQ_LOG] = res[IBAQ_NORMALIZED].apply(lambda x: (math.log10(x) + 10) if x > 0 else 0)

    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res[IBAQ_PPB] = res[IBAQ_NORMALIZED].apply(lambda x: x * 100000000)

    return res


def parse_uniprot_name(identifier: str) -> str:
    """
    Parse the uniprot name from the identifier  (e.g. sp|P12345|PROT_NAME)
    :param identifier: Uniprot identifier
    :return:
    """
    return identifier.split('|')[2]


@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option("-p", "--peptides",
              help="Peptide identifications with intensities following the peptide intensity output")
@click.option("-e", "--enzyme", help="Enzyme used during the analysis of the dataset (default: Trypsin)",
              default="Trypsin")
@click.option("-n", "--normalize", help="Normalize IBAQ values using by using the total IBAQ of the experiment",
              is_flag=True)
@click.option("--min_aa", help="Minimum number of amino acids to consider a peptide", default=7)
@click.option("--max_aa", help="Maximum number of amino acids to consider a peptide", default=30)
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
def ibaq_compute(fasta: str, peptides: str, enzyme: str, normalize: bool, min_aa: int, max_aa: int,
                 output: str) -> None:
    """
    This command computes the IBAQ values for a file output of peptides with the format described in
    peptide_contaminants_file_generation.py.
    :param min_aa: Minimum number of amino acids to consider a peptide
    :param max_aa: Maximum number of amino acids to consider a peptide
    :param fasta: Fasta file used to perform the peptide identification
    :param peptides: Peptide intensity file
    :param enzyme: Enzyme used to digest the protein sample
    :param normalize: use some basic normalization steps.
    :param output: output format containing the ibaq values.
    :return:
    """
    if peptides is None or fasta is None:
        print_help_msg(ibaq_compute)
        exit(1)

    fasta_proteins = list()  # type: list[FASTAEntry]
    FASTAFile().load(fasta, fasta_proteins)
    uniquepepcounts = dict()  # type: dict[str, int]
    digestor = ProteaseDigestion()
    digestor.setEnzyme(enzyme)

    def get_average_nr_peptides_unique_bygroup(pdrow: Series) -> Series:
        """
        Get the average intensity for protein groups
        :param pdrow: peptide row
        :return: average intensity
        """
        proteins = pdrow.name[0].split(';')
        summ = 0
        for prot in proteins:
            summ += uniquepepcounts[prot]
        if len(proteins) > 0 and summ > 0:
            return pdrow.NormIntensity / (summ / len(proteins))
        # If there is no protein in the group, return np nan
        return np.nan  # type: ignore

    for entry in fasta_proteins:
        digest = list()  # type: list[str]
        digestor.digest(AASequence().fromString(entry.sequence), digest, min_aa, max_aa)
        digestuniq = set(digest)
        uniquepepcounts[parse_uniprot_name(entry.identifier)] = len(digestuniq)

    data = pd.read_csv(peptides, sep=",")
    print(data.head())
    # next line assumes unique peptides only (at least per indistinguishable group)

    res = pd.DataFrame(data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()).apply(
        get_average_nr_peptides_unique_bygroup, 1)
    res = res.sort_values(ascending=False)
    res = res.to_frame()
    res = res.reset_index()
    res = res.rename(columns={0: IBAQ})

    if normalize:
        res = normalize_ibaq(res)
        # Remove IBAQ_NORMALIZED NAN values
        res = res.dropna(subset=[IBAQ_NORMALIZED])
        plot_column = IBAQ_PPB
    else:
        # Remove IBAQ NAN values
        res = res.dropna(subset=[IBAQ])
        plot_column = IBAQ

    plot_distributions(res, plot_column, SAMPLE_ID, log2=True)
    plot_box_plot(res, plot_column, SAMPLE_ID, log2=True,
                  title="IBAQ Distribution", violin=False)

    # # For absolute expression the relation is one sample + one condition
    # condition = data[CONDITION].unique()[0]
    # res[CONDITION] = condition.lower()

    res.to_csv(output, index=False)


if __name__ == '__main__':
    ibaq_compute()
