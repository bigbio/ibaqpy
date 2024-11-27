from unittest import TestCase
from ibaqpy.bin.compute_tpa import compute_tpa
from .common import datafile

class TestIbaqpy(TestCase):

    def test_ibaq_compute(self):
        args = {
            "fasta": datafile("Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta"),
            "peptides": datafile("PXD017834-peptides.csv"),
            "organism": "human",
            "ruler": True,
            "ploidy": 2,
            "cpc": 200,
            "output": datafile("PXD017834-tpa-norm.csv"),
            "verbose": True,
            "qc_report": datafile("TPA-QCprofile.pdf"),
        }
        print(args)
        compute_tpa(**args)
