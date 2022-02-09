# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import click
import pandas as pd
from pyopenms import *


def print_help_msg(command):
  with click.Context(command) as ctx:
    click.echo(command.get_help(ctx))

@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option("-p", "--peptides", help="Peptide identifications with intensities following the triqler output")
@click.option("-e", "--enzyme", help="Enzyme used during the analysis of the dataset (default: Trypsin)", default="Trypsin")
def ibaq_compute( fasta, peptides, enzyme):
  if peptides is None or fasta is None:
    print_help_msg(ibaq_compute)
    exit(-1)

  fasta_proteins = list()  # type: list[FASTAEntry]
  FASTAFile().load(fasta, fasta_proteins)
  uniquepepcounts = dict()  # type: dict[str, int]
  MINLEN = 6
  MAXLEN = 30
  ENZYMENAME = enzyme
  digestor = ProteaseDigestion()
  digestor.setEnzyme(ENZYMENAME)

  def get_average_nr_peptides_unique_bygroup(pdrow):
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
  ## next line assumes unique peptides only (at least per indistinguishable group)
  res = pd.DataFrame(data.groupby('proteins')['intensity'].sum()).apply(get_average_nr_peptides_unique_bygroup, 1)
  res.to_csv("res.csv")

if __name__ == '__main__':

  ibaq_compute()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
