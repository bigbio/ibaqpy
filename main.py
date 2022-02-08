# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd
from pyopenms import *

if __name__ == '__main__':
    fasta = list()  # type: list[FASTAEntry]
    FASTAFile().load("/Volumes/Data/Homo-sapiens-uniprot-reviewed-isoforms-contaminants-decoy-202105.fasta", fasta)
    uniquepepcounts = dict()  # type: dict[str, int]
    MINLEN = 6
    MAXLEN = 30
    ENZYMENAME = "Trypsin"
    digestor = ProteaseDigestion()
    digestor.setEnzyme(ENZYMENAME)


    def getAverageNrUniqPepsForIndistGrp(pdrow):
        proteins = pdrow.name.split(';')
        summ = 0
        for prot in proteins:
            summ += uniquepepcounts[prot]

        return pdrow.intesity / (summ / len(proteins))

    for entry in fasta:
        digest = list()  # type: list[str]
        digestor.digest(AASequence().fromString(entry.sequence), digest, MINLEN, MAXLEN)
        digestuniq = set(digest)
        uniquepepcounts[entry.identifier] = len(digestuniq)

    data = pd.read_csv("/Users/pfeuffer/Downloads/out_triqler_Sample1.tsv", sep="\t")
    print(data.head())
    ## next line assumes unique peptides only (at least per indistinguishable group)
    res = pd.DataFrame(data.groupby('proteins')['intensity'].sum()).apply(getAverageNrUniqPepsForIndistGrp, 1)
    res.to_csv("res.csv")


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
