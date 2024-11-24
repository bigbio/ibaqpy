from unittest import TestCase
from ibaqpy.ibaq.peptide_normalization import peptide_normalization
from ibaqpy.ibaq.compute_ibaq import ibaq_compute
from .common import datafile

class TestIbaqpy(TestCase):
    def test_feature_assembly(self):
        args = {
            "parquet": datafile("feature.parquet"),
            "sdrf": datafile("PXD017834-TMT.sdrf.tsv"),
            "min_aa": 7,
            "min_unique": 2,
            "remove_ids": None,
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": True,
            "output": datafile("PXD017834-peptides-norm.csv"),
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "max_min",
            "log2": True,
            "save_parquet": True,
        }
        print(__package__)
        print(args)
        peptide_normalization(**args)

    def test_ibaq_compute(self):
        args = {
            "fasta": datafile("Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta"),
            "peptides": datafile("PXD017834-peptides-norm.csv"),
            "enzyme": "Trypsin",
            "normalize": True,
            "min_aa": 7,
            "max_aa": 30,
            "output": datafile("PXD017834-ibaq-norm.csv"),
            "verbose": True,
            "qc_report": datafile("IBAQ-QCprofile.pdf"),
        }
        print(args)
        ibaq_compute(**args)
