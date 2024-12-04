from ibaqpy.bin.peptides2protein import peptides_to_protein

def test_ibaq_compute():
    args = {
        "fasta": "tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
        "peptides": "tests/example/PXD017834-peptides.csv",
        "enzyme": "Trypsin",
        "normalize": True,
        "min_aa": 7,
        "max_aa": 30,
        "tpa": True,
        "ruler": True,
        "ploidy": 2,
        "cpc": 200,
        "organism": "human",
        "output": "tests/example/PXD017834-norm.csv",
        "verbose": True,
        "qc_report": "tests/example/QCprofile.pdf",
    }
    print(args)
    peptides_to_protein(**args)
