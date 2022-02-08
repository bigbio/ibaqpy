from pyopenms import *

if __name__ == '__main__':
    fasta = list()  # type: list[FASTAEntry]
    FASTAFile().load("/Volumes/Data/Homo-sapiens-uniprot-reviewed-isoforms-contaminants-decoy-202105.fasta", fasta)
    uniquepepcounts = dict()  # type: dict[str, int]
    MINLEN = 6
    MAXLEN = 30
    ENZYMENAME = "Trypsin"