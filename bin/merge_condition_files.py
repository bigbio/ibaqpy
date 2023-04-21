#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os

import click
import pandas as pd

from ibaq.ibaqpy_commons import *


def print_help_msg(command) -> None:
    """
  Print help information
  :param command: command to print helps
  :return: print help
  """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


@click.command()
@click.option("-i", "--input", help="Folder with all the Intensity files", required=True)
@click.option("-o", "--output", help="Prefix name for the file to group by condition")
@click.option("-p", "--pattern", help="Prefix of the pattern name for all the files in the folder")
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
        g.to_csv(f'{output}/{k}-grouped-Intensities.csv', index=False)  # '{}.csv'.format(k)


if __name__ == '__main__':
    merge_condition_generation()
