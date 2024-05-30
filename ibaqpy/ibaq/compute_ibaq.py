import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame, Series
from pyopenms import *

from ibaqpy.ibaq.ibaqpy_commons import (
    CONDITION,
    IBAQ,
    IBAQ_LOG,
    IBAQ_NORMALIZED,
    IBAQ_PPB,
    NORM_INTENSITY,
    PROTEIN_NAME,
    SAMPLE_ID,
    plot_box_plot,
    plot_distributions,
    get_accession,
)


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
    res[IBAQ_LOG] = res[IBAQ_NORMALIZED].apply(
        lambda x: (math.log10(x) + 10) if x > 0 else 0
    )

    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res[IBAQ_PPB] = res[IBAQ_NORMALIZED].apply(lambda x: x * 100000000)

    return res


def ibaq_compute(
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
    if peptides is None or fasta is None:
        raise ValueError("Fasta and peptides files are required")

    fasta_proteins = list()  # type: list[FASTAEntry]
    protein_accessions = list()
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
        nonlocal map_size
        proteins = pdrow.name[0].split(";")
        summ = 0
        for prot in proteins:
            summ += uniquepepcounts[prot]
        if len(proteins) > 0 and summ > 0:
            return pdrow.NormIntensity / map_size[pdrow.name] / (summ / len(proteins))
        # If there is no protein in the group, return np nan
        return np.nan  # type: ignore

    for entry in fasta_proteins:
        digest = list()  # type: list[str]
        digestor.digest(AASequence().fromString(entry.sequence), digest, min_aa, max_aa)
        digestuniq = set(digest)
        # TODO: Try to get protein accessions from multiple databases.
        protein_name = get_accession(entry.identifier)
        uniquepepcounts[protein_name] = len(digestuniq)
        protein_accessions.append(protein_name)

    data = pd.read_csv(peptides, sep=",")
    data = data[data[PROTEIN_NAME].isin(protein_accessions)]
    print(data.head())
    # next line assumes unique peptides only (at least per indistinguishable group)
    map_size = data.groupby(["ProteinName", "SampleID", "Condition"]).size().to_dict()
    res = pd.DataFrame(
        data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()
    ).apply(get_average_nr_peptides_unique_bygroup, 1)
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

    # Print the distribution of the protein IBAQ values
    if verbose:
        plot_width = len(set(res["SampleID"])) * 0.5 + 10
        pdf = PdfPages(qc_report)
        density = plot_distributions(
            res,
            plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="IBAQ Distribution",
        )
        plt.show()
        pdf.savefig(density)
        box = plot_box_plot(
            res,
            plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="IBAQ Distribution",
            violin=False,
        )
        plt.show()
        pdf.savefig(box)
        pdf.close()

    # # For absolute expression the relation is one sample + one condition
    # condition = data[CONDITION].unique()[0]
    # res[CONDITION] = condition.lower()

    res.to_csv(output, index=False)
