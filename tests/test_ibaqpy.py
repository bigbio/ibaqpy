from unittest import TestCase
from ibaqpy.bin.peptides2protein import peptides_to_protein
from .common import datafile

class TestIbaqpy(TestCase):

    def test_ibaq_compute(self):
        args = {
            "fasta": datafile("Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta"),
            "peptides": datafile("PXD017834-peptides.csv"),
            "enzyme": "Trypsin",
            "normalize": True,
            "min_aa": 7,
            "max_aa": 30,
            "tpa": True,
            "ruler": True,
            "ploidy": 2,
            "cpc": 200,
            "organism": "human",
            "output": datafile("PXD017834-norm.csv"),
            "verbose": True,
            "qc_report": datafile("QCprofile.pdf"),
        }
        print(args)
        peptides_to_protein(**args)
