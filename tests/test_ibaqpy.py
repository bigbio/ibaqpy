from unittest import TestCase

from ibaqpy.ibaq.peptide_normalization import peptide_normalization
from ibaqpy.ibaq.compute_ibaq import ibaq_compute


class TestIbaqpy(TestCase):
    def test_feature_assembly(self):
        args = {
            "parquet": __package__ + "PXD003947/PXD003947-featrue.parquet",
            "sdrf": __package__ + "PXD003947/PXD003947.sdrf.tsv",
            "min_aa": 7,
            "min_unique": 2,
            "remove_ids": __package__ + "../data/contaminants_ids.tsv",
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": True,
            "output": __package__ + "PXD003947/PXD003947-peptides-norm.csv",
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "max_min",
            "log2": True,
            "save_parquet": True,
        }
        print(args)
        peptide_normalization(**args)

    def test_ibaq_compute(self):
        args = {
            "fasta": __package__
            + "PXD003947/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
            "peptides": __package__ + "PXD003947/PXD003947-peptides-norm.csv",
            "enzyme": "Trypsin",
            "normalize": True,
            "min_aa": 7,
            "max_aa": 30,
            "output": __package__ + "PXD003947/PXD003947-ibaq-norm.csv",
            "verbose": True,
            "qc_report": __package__ + "PXD003947/IBAQ-QCprofile.pdf",
        }
        print(args)
        ibaq_compute(**args)
