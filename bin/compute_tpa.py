#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import pandas as pd
from pyopenms import *

from bin.compute_ibaq import print_help_msg, parse_uniprot_name
from ibaq.ibaqpy_commons import PROTEIN_NAME, NORM_INTENSITY, SAMPLE_ID, CONDITION
from ibaq.ibaqpy_commons import plot_distributions, plot_box_plot
import numpy as np
import os

@click.command()
@click.option("-f", "--fasta", help="Protein database")
@click.option("-p", "--peptides", help="Peptide identifications with intensities following the peptide intensity output")
@click.option("-r", "--ruler", help="Whether to use proteomicRuler", default=True)
@click.option("-n", "--ploidy", help="ploidy number", default=2)
@click.option("-c", "--cpc", help="cellular protein concentration(g/L)", default=200)
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
def tpa_compute(fasta: str, peptides: str, ruler: bool, ploidy: int, cpc: float, output: str) -> None:
    """
    This command computes the protein copies and concentrations according to a file output of peptides with the
    format described in peptide_contaminants_file_generation.py.
    :param fasta: Fasta file used to perform the peptide identification
    :param peptides: Peptide intensity file.
    :param ruler: Whether to compute protein copies, weight and concentration.
    :param ploidy: ploidy number.
    :param cpc: cellular protein concentration(g/L).
    :param output: output format containing the ibaq values.
    :return:
    """
    if peptides is None or fasta is None:
        print_help_msg(tpa_compute)
        exit(1)

    data = pd.read_csv(peptides, sep=",")
    print(data.head())

    res = pd.DataFrame(data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum())
    res = res.reset_index()
    proteins = res["ProteinName"].unique().tolist()
    proteins = sum([i.split(";") for i in proteins], [])

    # calculate molecular weight of quantified proteins
    mw_dict = dict()
    fasta_proteins = list()  # type: list[FASTAEntry]
    FASTAFile().load(fasta, fasta_proteins)
    for entry in fasta_proteins:
        accession, name = entry.identifier.split("|")[1:]
        if name in proteins:
            mw = AASequence().fromString(entry.sequence).getMonoWeight()
            mw_dict[name] = mw

    # calculate TPA for every protein group
    def get_protein_group_mw(group: str) -> float:
        mw_list = [mw_dict[i] for i in group.split(";")]
        return sum(mw_list)

    res["MolecularWeight"] = res.apply(lambda x: get_protein_group_mw(x["ProteinName"]), axis=1)
    res["MolecularWeight"] = res["MolecularWeight"].fillna(1)
    res["MolecularWeight"] = res["MolecularWeight"].replace(0, 1)
    res["TPA"] = res["NormIntensity"] / res["MolecularWeight"]
    plot_distributions(res, "TPA", SAMPLE_ID, log2=True)
    plot_box_plot(res, "TPA", SAMPLE_ID, log2=True, title="TPA Distribution", violin=False)

    # calculate protein weight(ng) and concentration(nM)
    if ruler:
        avogadro = 6.02214129e23
        average_base_pair_mass = 617.96  # 615.8771

        organism = res.loc[0, "ProteinName"].split("_")[1].lower()
        histone_df = pd.read_json(open(os.path.split(__file__)[0] + os.sep + "histones.json", "rb")).T
        target_histones = histone_df[histone_df["name"] == organism.lower()]
        genome_size = target_histones["genome_size"].values[0]
        histones_list = target_histones["histone_entries"].values[0]                                                                          
        dna_mass = ploidy * genome_size * average_base_pair_mass / avogadro

        def calculate(protein_intensity, histone_intensity, mw):
            copy = (protein_intensity / histone_intensity) * dna_mass * avogadro / mw
            # The number of moles is equal to the number of particles divided by Avogadro's constant
            moles = copy * 1e9 / avogadro  # unit nmol
            weight = moles * mw  # unit ng
            return tuple([copy, moles, weight])

        def proteomicRuler(df):
            histone_intensity = df[df["ProteinName"].isin(histones_list)]["NormIntensity"].sum()
            histone_intensity = histone_intensity if histone_intensity > 0 else 1
            df[["Copy", "Moles[nmol]", "Weight[ng]"]] = df.apply(lambda x: calculate(x["NormIntensity"], histone_intensity, x["MolecularWeight"]), axis = 1, result_type="expand")
            volume = df["Weight[ng]"].sum() * 1e-9 / cpc  # unit L
            df["Concentration[nM]"] = df["Moles[nmol]"] / volume  # unit nM
            return df

        res = res.groupby(["Condition"]).apply(proteomicRuler)

        plot_distributions(res, "Concentration[nM]", SAMPLE_ID, log2=True)
        plot_box_plot(res, "Concentration[nM]", SAMPLE_ID, log2=True, title="Concentration[nM] Distribution", violin=False)
        res.to_csv(output, index=False)


if __name__ == '__main__':
    tpa_compute()
