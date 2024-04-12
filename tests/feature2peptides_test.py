from unittest import TestCase

from commands.features2peptides import peptide_normalization


class TestFeatureAssembly(TestCase):
    def test__load_sdrf_info(self):
        args = {
            "parquet": __package__ + "tests/PXD003947-featrue.parquet",
            "sdrf": __package__ + "tests/PXD003947.sdrf.tsv",
            "min_aa": 7,
            "min_unique": 2,
            "remove_ids": __package__ + "data/contaminants_ids.tsv",
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": True,
            "output": __package__ + "tests/PXD003947-peptides-norm.csv",
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "max_min",
            "log2": True,
            "save_parquet": True,
        }
        peptide_normalization(**args)
        
