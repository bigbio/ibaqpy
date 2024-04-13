#!/usr/bin/env python
# -*- coding: utf-8 -*-
import click
from ibaq.compute_ibaq import ibaq_compute


@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option(
    "-p",
    "--peptides",
    help="Peptide identifications with intensities following the peptide intensity output",
)
@click.option(
    "-e",
    "--enzyme",
    help="Enzyme used during the analysis of the dataset (default: Trypsin)",
    default="Trypsin",
)
@click.option(
    "-n",
    "--normalize",
    help="Normalize IBAQ values using by using the total IBAQ of the experiment",
    is_flag=True,
)
@click.option(
    "--min_aa", help="Minimum number of amino acids to consider a peptide", default=7
)
@click.option(
    "--max_aa", help="Maximum number of amino acids to consider a peptide", default=30
)
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
@click.option(
    "--verbose",
    help="Print addition information about the distributions of the intensities, number of peptides remove "
    "after normalization, etc.",
    is_flag=True,
)
@click.option(
    "--qc_report",
    help="PDF file to store multiple QC images",
    default="IBAQ-QCprofile.pdf",
)
def peptides2proteins(
    fasta: str,
    peptides: str,
    enzyme: str,
    normalize: bool,
    min_aa: int,
    max_aa: int,
    output: str,
    verbose: bool,
    qc_report: str,
) -> None:
    """
    This command computes the IBAQ values for a file output of peptides with the format described in
    peptide_contaminants_file_generation.py.
    :param min_aa: Minimum number of amino acids to consider a peptide.
    :param max_aa: Maximum number of amino acids to consider a peptide.
    :param fasta: Fasta file used to perform the peptide identification.
    :param peptides: Peptide intensity file.
    :param enzyme: Enzyme used to digest the protein sample.
    :param normalize: use some basic normalization steps.
    :param output: output format containing the ibaq values.
    :param verbose: Print addition information.
    :param qc_report: PDF file to store multiple QC images.
    :return:
    """
    ibaq_compute(**click.get_current_context().params)

if __name__ == "__main__":
    peptides2proteins()
