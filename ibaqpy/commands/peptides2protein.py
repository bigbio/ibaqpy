import click
from ibaqpy.ibaq.peptides2protein import peptides_to_protein
from ibaqpy.model.organism_metadata import OrganismDescription


@click.command("peptides2protein", short_help="Compute IBAQ values for proteins")
@click.option(
    "-f",
    "--fasta",
    help="Protein database to compute IBAQ values",
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
@click.option("--min_aa", help="Minimum number of amino acids to consider a peptide", default=7)
@click.option("--max_aa", help="Maximum number of amino acids to consider a peptide", default=30)
@click.option("-t", "--tpa", help="Whether calculate TPA", is_flag=True)
@click.option("-r", "--ruler", help="Whether to use ProteomicRuler", is_flag=True)
@click.option("-i", "--ploidy", help="Ploidy number (default: 2)", default=2)
@click.option(
    "-m",
    "--organism",
    help="Organism source of the data (default: human)",
    type=click.Choice(sorted(map(str.lower, OrganismDescription.registered_organisms())), case_sensitive=False),
    default="human")
@click.option("-c", "--cpc", help="Cellular protein concentration(g/L) (default: 200)", default=200)
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
    default="QCprofile.pdf",
)
@click.pass_context
def peptides2protein(
    click_context,
    fasta: str,
    peptides: str,
    enzyme: str,
    normalize: bool,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    ruler: bool,
    organism: str,
    ploidy: int,
    cpc: float,
    output: str,
    verbose: bool,
    qc_report: str,
) -> None:
    """
    This command computes the IBAQ/TPA values for a file output of peptides with the format described in
    peptide_contaminants_file_generation.py.
    :param click_context: Click context
    :param min_aa: Minimum number of amino acids to consider a peptide.
    :param max_aa: Maximum number of amino acids to consider a peptide.
    :param fasta: Fasta file used to perform the peptide identification.
    :param peptides: Peptide intensity file.
    :param enzyme: Enzyme used to digest the protein sample.
    :param normalize: use some basic normalization steps.
    :param tpa: Whether calculate TPA.
    :param ruler: Whether to compute protein copies, weight and concentration.
    :param organism: Organism source of the data.
    :param ploidy: Ploidy number.
    :param cpc: Cellular protein concentration(g/L).
    :param output: output format containing the ibaq values.
    :param verbose: Print addition information.
    :param qc_report: PDF file to store multiple QC images.
    :return:
    """
    peptides_to_protein(
        fasta=fasta,
        peptides=peptides,
        enzyme=enzyme,
        normalize=normalize,
        min_aa=min_aa,
        max_aa=max_aa,
        tpa=tpa,
        ruler=ruler,
        ploidy=ploidy,
        cpc=cpc,
        organism=organism,
        output=output,
        verbose=verbose,
        qc_report=qc_report,
    )
