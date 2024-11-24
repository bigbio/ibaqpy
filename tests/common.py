from pathlib import Path

TEST_DATA_ROOT = Path(__file__).parent / "example"
DATA_ROOT = Path(__file__).parent.parent / "data"

def datafile(path: str):
    path = str(path).removeprefix("/").removeprefix("example/")
    print(TEST_DATA_ROOT)
    return str(TEST_DATA_ROOT / path)