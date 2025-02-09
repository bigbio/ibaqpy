import logging

from ibaqpy.ibaq.peptide_normalization import peptide_normalization
from pathlib import Path

TESTS_DIR = Path(__file__).parent

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def test_feature_assembly():
    """
    Test the peptide normalization process by setting up arguments for the
    `peptide_normalization` function and executing it. This test checks the
    function's ability to process a feature table from a parquet file and an
    SDRF file, applying various filtering and normalization steps, and saving
    the output to a CSV file. It ensures that the output file is removed before
    the test to avoid conflicts.

    The test uses the following parameters:
    - parquet: Path to the input parquet file containing feature data.
    - sdrf: Path to the SDRF file for experimental metadata.
    - min_aa: Minimum number of amino acids required for peptides.
    - min_unique: Minimum number of unique peptides required for proteins.
    - remove_ids: Path to a file with protein IDs to remove, if any.
    - remove_decoy_contaminants: Flag to remove decoy and contaminant proteins.
    - remove_low_frequency_peptides: Flag to remove low-frequency peptides.
    - output: Path to the output CSV file for normalized peptide intensities.
    - skip_normalization: Flag to skip the normalization process.
    - nmethod: Method for feature-level normalization.
    - pnmethod: Method for peptide-level normalization.
    - log2: Flag to apply log2 transformation to intensities.
    - save_parquet: Flag to save the output as a parquet file.
    """

    args = {
        "parquet": str(TESTS_DIR / "example/feature.parquet"),
        "sdrf": str(TESTS_DIR / "example/PXD017834-TMT.sdrf.tsv"),
        "min_aa": 7,
        "min_unique": 2,
        "remove_ids": None,
        "remove_decoy_contaminants": True,
        "remove_low_frequency_peptides": True,
        "output": str(TESTS_DIR / "example" / "out" / "PXD017834-peptides-norm.csv"),
        "skip_normalization": False,
        "nmethod": "median",
        "pnmethod": "none",
        "log2": True,
        "save_parquet": True,
    }
    logger.info(args)
    out = Path(args["output"])
    if out.exists():
        out.unlink()
    peptide_normalization(**args)
