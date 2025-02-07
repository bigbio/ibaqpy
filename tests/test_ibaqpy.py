from ibaqpy.ibaq.peptides2protein import peptides_to_protein

from pathlib import Path

TESTS_DIR = Path(__file__).parent


def test_ibaq_compute():
    """
    Test the computation of IBAQ values using the peptides_to_protein function.

    This test sets up the necessary arguments, including paths to input files,
    enzyme type, normalization options, and output paths, and then calls the
    peptides_to_protein function to perform the computation. It prints the
    arguments for verification before execution.

    The test uses example data files located in the 'example' directory and
    outputs results to the 'out' directory.
    """
    args = {
        "fasta": str(
            TESTS_DIR / "example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta"
        ),
        "peptides": str(TESTS_DIR / "example/PXD017834-peptides.csv"),
        "enzyme": "Trypsin",
        "normalize": True,
        "min_aa": 7,
        "max_aa": 30,
        "tpa": True,
        "ruler": True,
        "ploidy": 2,
        "cpc": 200,
        "organism": "human",
        "output": str(TESTS_DIR / "example" / "out" / "PXD017834-ibaq.tsv"),
        "verbose": True,
        "qc_report": str(TESTS_DIR / "example/out/QCprofile.pdf"),
    }
    print(args)
    peptides_to_protein(**args)
