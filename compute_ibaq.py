import math

import click
import pandas as pd
from pandas import DataFrame, Series
from pyopenms import *

from ibaqpy_commons import remove_contaminants_decoys


def print_help_msg(command):
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


def normalize_ibaq(res: DataFrame) -> DataFrame:
    """
    Normalize the ibaq values using the total ibaq of the sample. The resulted
    ibaq values are then multiplied by 100'000'000 (PRIDE database noramalization)
    for the ibaq ppb and log10 shifted by 10 (ProteomicsDB)
    :param res: Dataframe
    :return:
    """

    # Get the total intensity for the sample.
    total_ibaq = res['ibaq'].sum()
    # Normalization method used by Proteomics DB 10 + log10(ibaq/sum(ibaq))
    res['ibaq_log'] = res['ibaq'].apply(lambda x: 100 + math.log2(x / total_ibaq))
    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res['ibaq_ppb'] = res['ibaq'].apply(lambda x: (x / total_ibaq) * 100000000)
    return res


@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option("-p", "--peptides", help="Peptide identifications with intensities following the triqler output")
@click.option("-e", "--enzyme", help="Enzyme used during the analysis of the dataset (default: Trypsin)",
              default="Trypsin")
@click.option("-n", "--normalize", help="Normalize IBAQ values using by using the total IBAQ of the experiment",
              is_flag=True)
@click.option("--contaminants_file", help="Contaminants protein accession", default="contaminants_ids.tsv")
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
def ibaq_compute(fasta: str, peptides: str, enzyme: str, normalize: bool, contaminants_file: str, output: str) -> None:
    if peptides is None or fasta is None:
        print_help_msg(ibaq_compute)
        exit(1)

    fasta_proteins = list()  # type: list[FASTAEntry]
    FASTAFile().load(fasta, fasta_proteins)
    uniquepepcounts = dict()  # type: dict[str, int]
    MINLEN = 6
    MAXLEN = 30
    ENZYMENAME = enzyme
    digestor = ProteaseDigestion()
    digestor.setEnzyme(ENZYMENAME)

    def get_average_nr_peptides_unique_bygroup(pdrow: Series) -> Series:
        proteins = pdrow.name.split(';')
        summ = 0
        for prot in proteins:
            summ += uniquepepcounts[prot]
        return pdrow.intensity / (summ / len(proteins))

    for entry in fasta_proteins:
        digest = list()  # type: list[str]
        digestor.digest(AASequence().fromString(entry.sequence), digest, MINLEN, MAXLEN)
        digestuniq = set(digest)
        uniquepepcounts[entry.identifier.decode('utf-8')] = len(digestuniq)

    data = pd.read_csv(peptides, sep="\t")
    print(data.head())
    # next line assumes unique peptides only (at least per indistinguishable group)

    res = pd.DataFrame(data.groupby('proteins')['intensity'].sum()).apply(get_average_nr_peptides_unique_bygroup, 1)
    res = res.sort_values(ascending=False)
    res = res.to_frame()
    res['protein'] = res.index
    res = res.rename(columns={0: "ibaq"})
    res = res[['protein', 'ibaq']]

    res = remove_contaminants_decoys(res, contaminants_file, protein_field='proteins')
    if normalize:
        res = normalize_ibaq(res)

    # For absolute expression the relation is one sample + one condition
    condition = data['condition'].unique()[0]
    res['condition'] = condition.lower()

    res.to_csv(output, index=False)


if __name__ == '__main__':
    ibaq_compute()
