import math

from typing import List, Union, Optional

import pandas as pd
import numpy as np

from pandas import DataFrame, Series
from matplotlib.backends.backend_pdf import PdfPages

from pyopenms import AASequence, ProteaseDigestion, FASTAFile

from ibaqpy.model.organism_metadata import OrganismDescription


from ibaqpy.ibaq.ibaqpy_commons import (
    CONDITION,
    IBAQ,
    IBAQ_LOG,
    IBAQ_NORMALIZED,
    IBAQ_PPB,
    NORM_INTENSITY,
    PROTEIN_NAME,
    SAMPLE_ID,
    TPA,
    MOLECULARWEIGHT,
    COPYNUMBER,
    CONCENTRATION_NM,
    MOLES_NMOL,
    WEIGHT_NG,
    is_parquet,
    plot_box_plot,
    plot_distributions,
    get_accession,
)


# Proteomic Ruler constants
AVAGADRO: float = 6.02214129e23
AVERAGE_BASE_PAIR_MASS: float = 617.96  # 615.8771


def normalize(group):
    """
    Normalize the ibaq values using the total ibaq of the sample.
    This method is called rIBAQ, originally published in https://pubs.acs.org/doi/10.1021/pr401017h
    :param group: Dataframe with all the ibaq values
    :return: Dataframe with the normalized ibaq values
    """
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
    res[IBAQ_PPB] = res[IBAQ_NORMALIZED].apply(lambda x: x * 100_000_000)
    return res


def handle_nonstandard_aa(aa_seq: str):
    """
    Any nonstandard amino acid will be removed.
    :param aa_seq: Protein sequences from multiple database.
    :return: One list contains nonstandard amoni acids and one remain sequence.
    """
    standard_aa = "ARNDBCEQZGHILKMFPSTWYV"
    nonstandard_aa_lst = [aa for aa in aa_seq if aa not in standard_aa]
    considered_seq = "".join([aa for aa in aa_seq if aa in standard_aa])
    return nonstandard_aa_lst, considered_seq


def extract_fasta(fasta: str, enzyme: str, proteins: List, min_aa: int, max_aa: int, tpa: bool):
    mw_dict = dict()
    fasta_proteins = list()
    FASTAFile().load(fasta, fasta_proteins)
    found_proteins = set()
    uniquepepcounts = dict()
    digestor = ProteaseDigestion()
    digestor.setEnzyme(enzyme)
    mw_dict = dict()
    for entry in fasta_proteins:
        accession = get_accession(entry.identifier)
        if accession in proteins:
            found_proteins.add(accession)
            digest = list()
            digestor.digest(AASequence().fromString(entry.sequence), digest, min_aa, max_aa)
            digestuniq = set(digest)
            uniquepepcounts[accession] = len(digestuniq)
            if tpa:
                try:
                    mw = AASequence().fromString(entry.sequence).getMonoWeight()
                    mw_dict[accession] = mw
                except ValueError:
                    error_aa, seq = handle_nonstandard_aa(entry.sequence)
                    mw = AASequence().fromString(seq).getMonoWeight()
                    mw_dict[accession] = mw
                    print(f"Nonstandard amino acids found in {accession}: {error_aa}, ignored!")
    if not found_proteins:
        raise ValueError(f"None of the {len(proteins)} proteins were found in the FASTA file")
    return uniquepepcounts, mw_dict, found_proteins


class ConcentrationWeightByProteomicRuler:
    """
    TODO

    Calculate protein copy number, moles, weight, and concentration for a given dataset
    for a given organism.

    This uses a proteomic ruler approach to estimate the copy number, moles,
    and weight of proteins in a dataset based on their normalized intensity and molecular
    weight. It also calculates the concentration in nM using the total weight and a
    provided concentration per cell (cpc).
    """
    organism: OrganismDescription
    ploidy: int
    concentration_per_cell: float
    dna_mass: float

    def __init__(self, organism: OrganismDescription, ploidy: int, concentration_per_cell: float):
        self.organism = organism
        self.ploidy = ploidy
        self.concentration_per_cell = concentration_per_cell

        self.dna_mass = (
            self.ploidy * self.organism.genome_size * AVERAGE_BASE_PAIR_MASS / AVAGADRO
        )

    def total_histone_intensities(self, protein_intensities: pd.DataFrame) -> float:
        histones = set(self.organism.histone_entries)
        is_histone_mask = protein_intensities[PROTEIN_NAME].isin(histones)
        histone_intensities = max(
            protein_intensities[is_histone_mask][NORM_INTENSITY].sum(),
            1.0
        )
        return histone_intensities

    def apply_ruler(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        histone_intensity = self.total_histone_intensities(protein_intensities)

        protein_intensities[COPYNUMBER] = (
            protein_intensities[NORM_INTENSITY]
            / histone_intensity
            * self.dna_mass
            * AVAGADRO
            / protein_intensities[MOLECULARWEIGHT]
        )

        protein_intensities[MOLES_NMOL] = protein_intensities[COPYNUMBER] * (1e9 / AVAGADRO)
        protein_intensities[WEIGHT_NG] = protein_intensities[MOLES_NMOL] * protein_intensities[MOLECULARWEIGHT]

        volume = protein_intensities[WEIGHT_NG].sum() / 1e-9 / self.concentration_per_cell
        protein_intensities[CONCENTRATION_NM] = volume * protein_intensities[MOLES_NMOL]
        return protein_intensities

    def __call__(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        return self.apply_ruler(protein_intensities)

    def apply_by_condition(self, protein_intensities: pd.DataFrame):
        protein_intensities = protein_intensities.groupby([CONDITION]).apply(self)
        return protein_intensities


class PeptideProteinMapper:
    _peptide_protein_ratio: dict[str, float]

    unique_peptide_counts: dict[str, int]
    map_size: dict[str, int]

    protein_mass_map: dict[str, float]

    def __init__(
        self,
        unique_peptide_counts: Optional[dict[str, int]] = None,
        map_size: Optional[dict[str, int]] = None,
        protein_mass_map: Optional[dict[str, float]] = None,
    ):
        self.unique_peptide_counts = unique_peptide_counts or {}
        self.map_size = map_size or {}
        self.protein_mass_map = protein_mass_map or {}

        self._peptide_protein_ratio = {}

    def peptide_protein_ratio(self, protein_group: str):
        if protein_group in self._peptide_protein_ratio:
            return self._peptide_protein_ratio[protein_group]

        proteins_list = protein_group.split(";")

        total = 0
        for prot in proteins_list:
            total += self.unique_peptide_counts[prot]

        if not proteins_list:
            val = self._peptide_protein_ratio[protein_group] = 0
        else:
            val = self._peptide_protein_ratio[protein_group] = total / len(proteins_list)
        return val

    def get_average_nr_peptides_unique_by_group(self, pdrow: Series) -> Union[float, Series]:
        """
        Calculate the average number of unique peptides per protein group.

        This function computes the average number of unique peptides for a given
        protein group by dividing the normalized intensity by the product of the
        group size and the average unique peptide count. If no proteins are found
        in the group or the sum of unique peptide counts is zero, it returns NaN.

        :param pdrow: Series containing the protein group data.
        :return: Series containing the average number of unique peptides per protein group.
        """

        average_peptides_per_protein = self.peptide_protein_ratio(pdrow.name[0])

        if average_peptides_per_protein > 0:
            return (
                pdrow.NormIntensity
                / self.map_size[pdrow.name]
                / average_peptides_per_protein
            )

        # If there is no protein in the group, return np nan
        return np.nan  # type: ignore

    def protein_group_mass(self, protein_group: str):
        """
        Calculate the molecular weight of a protein group.

        :param group: Protein group.
        :return: Molecular weight of the protein group.
        """
        mw_list = [self.protein_mass_map[i] for i in protein_group.split(";")]
        return sum(mw_list)


def peptides_to_protein(
    fasta: str,
    peptides: str,
    enzyme: str,
    normalize: bool,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    ruler: bool,
    ploidy: int,
    cpc: float,
    organism: str,
    output: str,
    verbose: bool,
    qc_report: str,
) -> None:
    """
    Compute IBAQ values for peptides and generate a QC report.

    This function processes peptide intensity data to compute IBAQ values,
    optionally normalizes the data, and calculates protein metrics such as
    weight and concentration using a proteomic ruler approach. It also
    generates a QC report with distribution plots if verbose mode is enabled.

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

    if organism:
        organism_descr = OrganismDescription.get(organism)
        if organism_descr is None:
            raise KeyError(f"Could not resolve organism description for {organism}")
    else:
        organism_descr = None

    if ruler:
        if not ploidy or not cpc or not organism or not tpa:
            raise ValueError(
                "Arguments `ploidy`, `cpc`, `organism` and `tpa` are required for calculate protein weight(ng) and concentration(nM)"
            )

    # load data
    if is_parquet(peptides):
        data = pd.read_parquet(peptides)
    else:
        data = pd.read_csv(peptides)
    data[NORM_INTENSITY] = data[NORM_INTENSITY].astype(float)
    data = data.dropna(subset=[NORM_INTENSITY])
    data = data[data[NORM_INTENSITY] > 0]

    # get fasta info
    proteins = data[PROTEIN_NAME].unique().tolist()

    # TODO: This should check for duplicates?
    proteins = sum([i.split(";") for i in proteins], [])

    unique_peptide_counts, mw_dict, found_proteins = extract_fasta(
        fasta, enzyme, proteins, min_aa, max_aa, tpa
    )

    data = data[data[PROTEIN_NAME].isin(found_proteins)]

    # data processing
    print(data.head())
    map_size = data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION]).size().to_dict()
    res = pd.DataFrame(data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum())

    protein_mapper = PeptideProteinMapper(
        unique_peptide_counts=unique_peptide_counts,
        map_size=map_size,
        protein_mass_map=mw_dict,
    )

    # ibaq
    res[IBAQ] = res.apply(protein_mapper.get_average_nr_peptides_unique_by_group, 1)
    res = res.reset_index()

    # normalize ibaq
    if normalize:
        res = normalize_ibaq(res)
        # Remove IBAQ_NORMALIZED NAN values
        res = res.dropna(subset=[IBAQ_NORMALIZED])
        plot_column = IBAQ_PPB
    else:
        # Remove IBAQ NAN values
        res = res.dropna(subset=[IBAQ])
        plot_column = IBAQ

    res = res.reset_index(drop=True)

    # tpa
    if tpa:
        res[MOLECULARWEIGHT] = (
            res[PROTEIN_NAME]
            .apply(protein_mapper.protein_group_mass)
            .fillna(1.0)
            .replace(0.0, 1.0)
        )
        res[TPA] = res[NORM_INTENSITY] / res[MOLECULARWEIGHT]

    # calculate protein weight(ng) and concentration(nM)
    if ruler:
        concentration_by_ruler = ConcentrationWeightByProteomicRuler(organism_descr, ploidy, cpc)
        res = concentration_by_ruler.apply_by_condition(res)

    # Print the distribution of the protein IBAQ values
    if verbose:
        plot_width = len(set(res[SAMPLE_ID])) * 0.5 + 10
        pdf = PdfPages(qc_report)
        density1 = plot_distributions(
            res,
            plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="{} Distribution".format(plot_column),
        )
        box1 = plot_box_plot(
            res,
            plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="{} Distribution".format(plot_column),
            violin=False,
        )
        pdf.savefig(density1, bbox_inches="tight")
        pdf.savefig(box1, bbox_inches="tight")
        if tpa:
            density2 = plot_distributions(
                res, TPA, SAMPLE_ID, log2=True, width=plot_width, title="TPA Distribution"
            )
            box2 = plot_box_plot(
                res,
                TPA,
                SAMPLE_ID,
                log2=True,
                width=plot_width,
                title="{} Distribution".format(TPA),
                violin=False,
            )
            pdf.savefig(density2, bbox_inches="tight")
            pdf.savefig(box2, bbox_inches="tight")
        if ruler:
            density3 = plot_distributions(
                res,
                COPYNUMBER,
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="{} Distribution".format(COPYNUMBER),
            )
            box3 = plot_box_plot(
                res,
                COPYNUMBER,
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="{} Distribution".format(COPYNUMBER),
                violin=False,
            )
            pdf.savefig(density3, bbox_inches="tight")
            pdf.savefig(box3, bbox_inches="tight")
            density4 = plot_distributions(
                res,
                CONCENTRATION_NM,
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="{} Distribution".format(CONCENTRATION_NM),
            )
            box4 = plot_box_plot(
                res,
                CONCENTRATION_NM,
                SAMPLE_ID,
                width=plot_width,
                log2=True,
                title="{} Distribution".format(CONCENTRATION_NM),
                violin=False,
            )
            pdf.savefig(density4, bbox_inches='tight')
            pdf.savefig(box4, bbox_inches="tight")
        pdf.close()

    res.to_csv(output, sep="\t", index=False)
