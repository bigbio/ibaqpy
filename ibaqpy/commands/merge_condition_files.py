
import os
import click
import pandas as pd

from ibaqpy.bin.ibaqpy_commons import CONDITION


@click.command()
@click.option(
    "-i", "--input", help="Folder with all the Intensity files", required=True
)
@click.option("-o", "--output", help="Prefix name for the file to group by condition")
@click.option(
    "-p", "--pattern", help="Prefix of the pattern name for all the files in the folder"
)
def merge_condition_generation(input: str, output: str, pattern: str) -> None:
    """
    Merge all the files in a folder with the specific pattern
    :param input: Input folder containing all the peptide Intensity files
    :param output: Output file prefix with all the intensities
    :param pattern: pattern of the files with the corresponding file name prefix
    :return:
    """

    files = [f for f in os.listdir(input) if pattern in f]
    df_from_each_file = (pd.read_csv(input + "/" + f) for f in files)
    concatenated_df = pd.concat(df_from_each_file, ignore_index=True)
    concatenated_df[CONDITION] = concatenated_df[CONDITION].str.lower()
    print(concatenated_df.head())

    for k, g in concatenated_df.groupby([CONDITION]):
        g.to_csv(
            f"{output}/{k}-grouped-Intensities.csv", index=False
        )  # '{}.csv'.format(k)
