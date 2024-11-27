import click
from ibaqpy.bin.compute_tpa import compute_tpa

@click.command("computetpa", short_help="Compute TPA values.")
@click.option(
    "-f",
    "--fasta",
    help="Protein database",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-p",
    "--peptides",
    help="Peptide identifications with intensities following the peptide intensity output",
    required=True,
    type=click.Path(exists=True),
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
    compute_tpa(
        fasta=fasta,
        peptides=peptides,
        ruler=ruler,
        organism=organism,
        ploidy=ploidy,
        cpc=cpc,
        output=output,
        verbose=verbose,
        qc_report=qc_report,
    )