from ibaqpy.bin.peptide_normalization import peptide_normalization


def test_feature_assembly():
    args = {
        "parquet": "tests/example/feature.parquet",
        "sdrf": "tests/example/PXD017834-TMT.sdrf.tsv",
        "min_aa": 7,
        "min_unique": 2,
        "remove_ids": None,
        "remove_decoy_contaminants": True,
        "remove_low_frequency_peptides": True,
        "output": "tests/example/PXD017834-peptides-norm.csv",
        "skip_normalization": False,
        "nmethod": "median",
        "pnmethod": "max_min",
        "log2": True,
        "save_parquet": True,
    }
    print(args)
    peptide_normalization(**args)