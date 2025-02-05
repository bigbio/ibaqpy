import logging
from pathlib import Path
import pytest

import pandas as pd

from ibaqpy.commands.correct_batches import run_batch_correction
TESTS_DIR = Path(__file__).parent


def test_correct_batches():
    # Test valid case
    args = {
        "folder": TESTS_DIR.parent / "data/ibaq-raw-hela",
        "pattern": "*ibaq.tsv",
        "comment": "#",
        "sep": "\t",
        "output": TESTS_DIR / "example/ibaq_corrected_combined.tsv",
        "sample_batch_map": None,
        "sample_id_column": "SampleID",
        "protein_id_column": "ProteinName",
        "ibaq_column": "Ibaq"
    }
    logging.debug("Arguments for run_batch_correction: %s", args)
    run_batch_correction(**args)
    
    # Assert the output file is created and not empty
    output_path = Path(args["output"])
    assert output_path.exists(), f"Expected output file {output_path} was not created."
    df = pd.read_csv(output_path, sep=args["sep"])
    assert not df.empty, "The corrected output file is empty."
    
    # Test invalid sample IDs
    with pytest.raises(ValueError):
        args["folder"] = TESTS_DIR.parent / "data/invalid-samples"
        run_batch_correction(**args)
    
    # Test missing required columns
    with pytest.raises(KeyError):
        args["folder"] = TESTS_DIR.parent / "data/ibaq-raw-hela"
        args["sample_id_column"] = "NonexistentColumn"
        run_batch_correction(**args)
    
    # Test invalid file pattern
    with pytest.raises(FileNotFoundError):
        args["pattern"] = "nonexistent*.tsv"
        run_batch_correction(**args)
if __name__ == "__main__":
    test_correct_batches()


