from ibaqpy.ibaq.peptide_normalization import peptide_normalization
from pathlib import Path

TESTS_DIR = Path(__file__).parent


def test_feature_assembly():
    args = {
        "parquet": str(TESTS_DIR / "example/feature.parquet"),
        "sdrf": str(TESTS_DIR / "example/PXD017834-TMT.sdrf.tsv"),
        "min_aa": 7,
        "min_unique": 2,
        "remove_ids": None,
        "remove_decoy_contaminants": True,
        "remove_low_frequency_peptides": True,
        "output": str(TESTS_DIR / "example/PXD017834-peptides-norm.csv"),
        "skip_normalization": False,
        "nmethod": "median",
        "pnmethod": "max_min",
        "log2": True,
        "save_parquet": True,
    }
    print(args)
    peptide_normalization(**args)
