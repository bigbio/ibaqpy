from unittest import TestCase
from ibaqpy.bin.peptide_normalization import peptide_normalization
from .common import datafile

class TestPeptide(TestCase):
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