import logging
import re
from pathlib import Path
from typing import Union

import click
import pandas as pd

from ibaqpy.ibaq.file_utils import create_anndata, combine_ibaq_tsv_files
from ibaqpy.ibaq.ibaqpy_commons import SAMPLE_ID_REGEX, SAMPLE_ID, PROTEIN_NAME, IBAQ, IBAQ_BEC
from ibaqpy.ibaq.ibaqpy_postprocessing import (
    pivot_wider,
    pivot_longer,
)
from ibaqpy.ibaq.utils import apply_batch_correction


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def is_valid_sample_id(
    samples: Union[str, list, pd.Series], sample_id_pattern: str = SAMPLE_ID_REGEX
) -> bool:
    """
    Validate sample IDs against a specified pattern.

    This function checks whether the provided sample IDs match a given regex pattern.
    It accepts a single sample ID, a list of sample IDs, or a pandas Series of sample IDs.
    If any sample ID does not match the pattern, it prints the invalid IDs and returns False.
    Otherwise, it returns True.

    Parameters
    ----------
    samples : Union[str, list, pd.Series]
        The sample ID(s) to validate.
    sample_id_pattern : str, optional
        The regex pattern to validate the sample IDs against. Defaults to 'SAMPLE_ID_REGEX'.

    Returns
    -------
    bool
        True if all sample IDs are valid, False otherwise.
    """
    sample_pattern = re.compile(sample_id_pattern)

    # Ensure samples is a list for uniform processing
    if isinstance(samples, str):
        samples = [samples]
    elif isinstance(samples, pd.Series):
        samples = samples.tolist()

    # Identify invalid sample names.
    invalid_samples = [sample for sample in samples if not sample_pattern.fullmatch(sample)]

    if invalid_samples:
        print("The following sample IDs are invalid:")
        for invalid_sample in invalid_samples:
            print(f" - {invalid_sample}")
        return False
    return True


def get_batch_id_from_sample_names(samples: list) -> list:
    """
    Extract batch IDs from a list of sample names.

    Each sample name is expected to have a batch ID as a prefix, separated by a hyphen.
    The function validates that the batch ID consists of alphanumeric characters only.
    Returns a list of unique batch IDs as integer factors.

    Parameters
    ----------
    samples : list
        A list of sample names, each containing a batch ID prefix.

    Returns
    -------
    list
        A list of integer factors representing unique batch IDs.

    Raises
    ------
    ValueError
        If a sample name does not contain a valid batch ID prefix or if the
        batch ID contains non-alphanumeric characters.
    """
    batch_ids = []
    for sample in samples:
        parts = sample.split("-")
        if not parts or not parts[0]:
            raise ValueError(f"Invalid sample name format: {sample}. Expected batch-id prefix.")
        batch_id = parts[0]
        if not re.match(r"^[A-Za-z0-9]+$", batch_id):
            raise ValueError(
                f"Invalid batch ID format: {batch_id}. Expected alphanumeric characters only."
            )
        batch_ids.append(batch_id)
    return pd.factorize(batch_ids)[0]


def run_batch_correction(
    folder: str,
    pattern: str,
    comment: str,
    sep: str,
    output: str,
    sample_id_column: str = SAMPLE_ID,
    protein_id_column: str = PROTEIN_NAME,
    ibaq_raw_column: str = IBAQ,
    ibaq_corrected_column: str = IBAQ_BEC,
    export_anndata: bool = False,
) -> pd.DataFrame:
    """
    Run batch correction on iBAQ data from TSV files in a specified directory.

    This function combines multiple TSV files, reshapes the data, validates sample IDs,
    applies batch correction, and optionally exports the results to an AnnData object.

    Parameters
    ----------
    folder : str
        Directory containing the TSV files.
    pattern : str
        Pattern to match files in the directory.
    comment : str
        Character indicating the start of a comment line in the TSV files.
    sep : str
        Delimiter for reading the TSV files.
    output : str
        File path to save the corrected iBAQ values.
    sample_id_column : str, optional
        Column name for sample IDs. Defaults to 'SAMPLE_ID'.
    protein_id_column : str, optional
        Column name for protein IDs. Defaults to 'PROTEIN_NAME'.
    ibaq_raw_column : str, optional
        Column name for raw iBAQ values. Defaults to 'IBAQ'.
    ibaq_corrected_column : str, optional
        Column name for corrected iBAQ values. Defaults to 'IBAQ_BEC'.
    export_anndata : bool, optional
        Whether to export the data to an AnnData object. Defaults to False.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the original and corrected iBAQ values.

    Raises
    ------
    ValueError
        If input files cannot be loaded, sample IDs are invalid, or output file cannot be saved.
    FileNotFoundError
        If the output file does not exist when exporting to AnnData.
    """

    # Load the data
    logger.info(f"Loading iBAQ data from TSV files in folder '{folder}'")

    try:
        df_ibaq = combine_ibaq_tsv_files(folder, pattern=pattern, comment=comment, sep=sep)
    except Exception as e:
        raise ValueError(f"Failed to load input files: {str(e)}")

    # Reshape the data to wide format
    df_wide = pivot_wider(
        df_ibaq,
        row_name=protein_id_column,
        col_name=sample_id_column,
        values=ibaq_raw_column,
        fillna=True,
    )

    # Validate the sample IDs
    if not is_valid_sample_id(df_wide.columns, SAMPLE_ID_REGEX):
        raise ValueError("Invalid sample IDs found in the data.")

    # Get the batch IDs
    batch_ids = get_batch_id_from_sample_names(df_wide.columns)

    # Run batch correction
    logger.info("Applying batch correction to iBAQ values")
    df_corrected = apply_batch_correction(df_wide, list(batch_ids), kwargs={})

    # Convert the data back to long format
    df_corrected = df_corrected.reset_index()
    df_corrected_long = pivot_longer(
        df_corrected,
        row_name=protein_id_column,
        col_name=sample_id_column,
        values=ibaq_corrected_column,
    )

    # Add the corrected ibaq values to the original dataframe.
    # Use sample/protein ID keys to merge the dataframes.
    df_ibaq = df_ibaq.merge(
        df_corrected_long, how="left", on=[sample_id_column, protein_id_column]
    )

    # Save the corrected iBAQ values to a file
    if output:
        try:
            df_ibaq.to_csv(output, sep=sep, index=False)
        except Exception as e:
            raise ValueError(f"Failed to save output file: {str(e)}")

    # Export the raw and corrected iBAQ values to an AnnData object
    if export_anndata:
        logger.info("Exporting raw and corrected iBAQ values to an AnnData object")
        output_path = Path(output)
        if not output_path.exists():
            raise FileNotFoundError(f"Output file {output} does not exist!")
        adata = create_anndata(
            df_ibaq,
            obs_col=sample_id_column,
            var_col=protein_id_column,
            value_col=ibaq_raw_column,
            layer_cols=[ibaq_corrected_column],
        )
        adata_filename = output_path.with_suffix(".h5ad")
        try:
            adata.write(adata_filename)
        except Exception as e:
            raise ValueError(f"Failed to write AnnData object: {e}")

    logger.info("Batch correction completed...")

    return df_ibaq


@click.command("correct-batches", short_help="Batch effect correction for iBAQ values.")
@click.option(
    "-f",
    "--folder",
    help="Folder that contains all TSV files with raw iBAQ values",
    required=True,
    default=None,
)
@click.option(
    "-p",
    "--pattern",
    help="Pattern for the TSV files with raw iBAQ values",
    required=True,
    default="*ibaq.tsv",
)
@click.option(
    "--comment",
    help="Comment character for the TSV files. Lines starting with this character will be ignored.",
    required=False,
    default="#",
)
@click.option("--sep", help="Separator for the TSV files", required=False, default="\t")
@click.option(
    "-o",
    "--output",
    help="Output file name for the combined iBAQ corrected values",
    required=True,
)
@click.option(
    "-sid",
    "--sample_id_column",
    help="Sample ID column name",
    required=False,
    default=SAMPLE_ID,
)
@click.option(
    "-pid",
    "--protein_id_column",
    help="Protein ID column name",
    required=False,
    default=PROTEIN_NAME,
)
@click.option(
    "-ibaq", "--ibaq_raw_column", help="Name of the raw iBAQ column", required=False, default=IBAQ
)
@click.option(
    "--ibaq_corrected_column",
    help="Name for the corrected iBAQ column",
    required=False,
    default=IBAQ_BEC,
)
@click.option(
    "--export_anndata",
    help="Export the raw and corrected iBAQ values to an AnnData object",
    is_flag=True,
)
@click.pass_context
def correct_batches(
    ctx,
    folder: str,
    pattern: str,
    comment: str,
    sep: str,
    output: str,
    sample_id_column: str,
    protein_id_column: str,
    ibaq_raw_column: str,
    ibaq_corrected_column: str,
    export_anndata: bool,
):
    """
    Correcting batch effects in iBAQ data.

    This command processes TSV files containing raw iBAQ values, applies batch correction,
    and outputs the corrected values. It supports various options for specifying file patterns,
    column names, and output formats, including exporting to an AnnData file.
    """
    run_batch_correction(
        folder=folder,
        pattern=pattern,
        comment=comment,
        sep=sep,
        output=output,
        sample_id_column=sample_id_column,
        protein_id_column=protein_id_column,
        ibaq_raw_column=ibaq_raw_column,
        ibaq_corrected_column=ibaq_corrected_column,
        export_anndata=export_anndata,
    )
