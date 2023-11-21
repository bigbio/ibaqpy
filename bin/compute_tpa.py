#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import click
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pyopenms import *

from bin.compute_ibaq import print_help_msg
from ibaq.ibaqpy_commons import (CONDITION, NORM_INTENSITY, PROTEIN_NAME, SAMPLE_ID,
                                 plot_box_plot, plot_distributions,
                                 remove_contaminants_decoys, get_accession)


def handle_nonstandard_aa(aa_seq: str) -> (list, str):
    """Any nonstandard amoni acid will be removed.

    :param aa_seq: Protein sequences from multiple database.
    :return: One list contains nonstandard amoni acids and one remain sequence.
    """
    standard_aa = 'ARNDBCEQZGHILKMFPSTWYV'
    nonstandard_aa_lst = [aa for aa in aa_seq if aa not in standard_aa]
    considered_seq = ''.join([aa for aa in aa_seq if aa in standard_aa])
    return nonstandard_aa_lst, considered_seq


@click.command()
@click.option("-f", "--fasta", help="Protein database")
@click.option(
    "-p",
    "--peptides",
    help="Peptide identifications with intensities following the peptide intensity output",
)
@click.option("-r", "--ruler", help="Whether to use ProteomicRuler", is_flag=True)
@click.option("-n", "--ploidy", help="Ploidy number", default=2)
@click.option("-m", "--organism", help="Organism source of the data", default="human")
@click.option("-c", "--cpc", help="Cellular protein concentration(g/L)", default=200)
@click.option("-o", "--output", help="Output file with the proteins and other values")
@click.option(
    "--verbose",
    help="Print addition information about the distributions of the intensities, number of peptides remove "
    "after normalization, etc.",
    is_flag=True,
)
@click.option(
    "--qc_report",
    help="PDF file to store multiple QC images",
    default="TPA-QCprofile.pdf",
)
def tpa_compute(
    fasta: str,
    peptides: str,
    ruler: bool,
    organism: str,
    ploidy: int,
    cpc: float,
    output: str,
    verbose: bool,
    qc_report: str,
) -> None:
    """
    This command computes the protein copies and concentrations according to a file output of peptides with the
    format described in peptide_contaminants_file_generation.py.
    :param fasta: Fasta file used to perform the peptide identification.
    :param peptides: Peptide intensity file without normalization.
    :param ruler: Whether to compute protein copies, weight and concentration.
    :param organism: Organism source of the data.
    :param ploidy: Ploidy number.
    :param cpc: Cellular protein concentration(g/L).
    :param output: Output format containing the TPA values, protein copy numbers and concentrations.
    :param verbose: Print addition information.
    :param qc_report: PDF file to store multiple QC images.
    :return:
    """
    if peptides is None or fasta is None:
        print_help_msg(tpa_compute)
        exit(1)

    data = pd.read_csv(
        peptides, sep=",", usecols=[PROTEIN_NAME, NORM_INTENSITY, SAMPLE_ID, CONDITION]
    )
    data[NORM_INTENSITY] = data[NORM_INTENSITY].astype("float")
    data = data.dropna(subset=[NORM_INTENSITY])
    data = data[data[NORM_INTENSITY] > 0]
    print(data.head())

    res = pd.DataFrame(
        data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()
    )
    res = res.reset_index()
    proteins = res[PROTEIN_NAME].unique().tolist()
    proteins = sum([i.split(";") for i in proteins], [])

    # calculate molecular weight of quantified proteins
    mw_dict = dict()
    fasta_proteins = list()  # type: list[FASTAEntry]
    FASTAFile().load(fasta, fasta_proteins)
    for entry in fasta_proteins:
        accession = get_accession(entry.identifier)
        if accession in proteins:
            try:
                mw = AASequence().fromString(entry.sequence).getMonoWeight()
                mw_dict[accession] = mw
            except:
                error_aa, seq = handle_nonstandard_aa(entry.sequence)
                mw = AASequence().fromString(seq).getMonoWeight()
                mw_dict[accession] = mw
                print(f"Nonstandard amimo acids found in {accession}: {error_aa}, ignored!")

    res = res[res[PROTEIN_NAME].isin(mw_dict.keys())]

    # calculate TPA for every protein group
    def get_protein_group_mw(group: str) -> float:
        mw_list = [mw_dict[i] for i in group.split(";")]
        return sum(mw_list)

    res["MolecularWeight"] = res.apply(
        lambda x: get_protein_group_mw(x[PROTEIN_NAME]), axis=1
    )
    res["MolecularWeight"] = res["MolecularWeight"].fillna(1)
    res["MolecularWeight"] = res["MolecularWeight"].replace(0, 1)
    res["TPA"] = res[NORM_INTENSITY] / res["MolecularWeight"]
    # Print the distribution of the protein TPA values
    if verbose:
        plot_width = len(set(res[SAMPLE_ID])) * 0.5 + 10
        pdf = PdfPages(qc_report)
        density = plot_distributions(
            res, "TPA", SAMPLE_ID, log2=True, width=plot_width, title="TPA Distribution"
        )
        plt.show()
        pdf.savefig(density)
        box = plot_box_plot(
            res,
            "TPA",
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="TPA Distribution",
            violin=False,
        )
        plt.show()
        pdf.savefig(box)

    # calculate protein weight(ng) and concentration(nM)
    if ruler:
        avogadro = 6.02214129e23
        average_base_pair_mass = 617.96  # 615.8771

        organism = organism.lower()
        histone_df = pd.read_json(
            open(os.path.split(__file__)[0] + os.sep + "histones.json", "rb")
        ).T
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

        def proteomic_ruler(df):
            histone_intensity = df[df[PROTEIN_NAME].isin(histones_list)][
                NORM_INTENSITY
            ].sum()
            histone_intensity = histone_intensity if histone_intensity > 0 else 1
            df[["Copy", "Moles[nmol]", "Weight[ng]"]] = df.apply(
                lambda x: calculate(
                    x[NORM_INTENSITY], histone_intensity, x["MolecularWeight"]
                ),
                axis=1,
                result_type="expand",
            )
            volume = df["Weight[ng]"].sum() * 1e-9 / cpc  # unit L
            df["Concentration[nM]"] = df["Moles[nmol]"] / volume  # unit nM
            return df

        res = res.groupby([CONDITION]).apply(proteomic_ruler)

        if verbose:
            density = plot_distributions(
                res, "Copy", SAMPLE_ID, width=plot_width, log2=True, title="Copy numbers Distribution"
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                res,
                "Copy",
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="Copy numbers Distribution",
                violin=False,
            )
            plt.show()
            pdf.savefig(box)

            density = plot_distributions(
                res,
                "Concentration[nM]",
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="Concentration[nM] Distribution",
            )
            plt.show()
            pdf.savefig(density)
            box = plot_box_plot(
                res,
                "Concentration[nM]",
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="Concentration[nM] Distribution",
                violin=False,
            )
            plt.show()
            pdf.savefig(box)
            pdf.close()
        res.to_csv(output, index=False)


if __name__ == "__main__":
    tpa_compute()
