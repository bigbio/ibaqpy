from pathlib import Path

TEST_DATA_ROOT = Path(__file__).parent / "example"

def datafile(path: str):
    path = str(path)
    return str(TEST_DATA_ROOT / path)