import logging
from pathlib import Path
import click
from ibaqpy.commands.features2peptides import features2parquet
from ibaqpy.commands.peptides2protein import peptides2protein
from ibaqpy.commands.tsne_visualization import tsne_visualization
from ibaqpy.commands.correct_batches import correct_batches

import ibaqpy.__init__ as __init__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


LOG_LEVELS = ["debug", "info", "warn"]
LOG_LEVELS_TO_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warn": logging.WARN,
}

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    version=__init__.__version__,
    package_name="ibaqpy",
    message="%(package)s %(version)s",
)
@click.option("-v", "--log-level", type=click.Choice(LOG_LEVELS, False), default="debug", help="Set the logging level.")
@click.option("--log-file", type=click.Path(writable=True, path_type=Path), required=False, help="Write log to this file.")
def cli(log_level: str, log_file: Path):
    """
    Aggregrate and normalize quantitative proteomics using iBAQ (Intensity-Based Absolute Quantification)
    for the quantms ecosystem.
    """

    logging.basicConfig(
        format="%(asctime)s [%(funcName)s] - %(message)s", level=LOG_LEVELS_TO_LEVELS[log_level.lower()]
    )
    logging.captureWarnings()

    if log_file:
        if not log_file.exists():
            if not log_file.parent.exists():
                log_file.parent.mkdir(parents=True, exist_ok=True)
            handler = logging.FileHandler(log_file)
            handler.setLevel(LOG_LEVELS_TO_LEVELS[log_level.lower()])
            handler.setFormatter(
                logging.Formatter("%(asctime)s [%(funcName)s] - %(message)s")
            )
            logging.getLogger().addHandler(handler)



cli.add_command(features2parquet)
cli.add_command(peptides2protein)
cli.add_command(tsne_visualization)
cli.add_command(correct_batches)


def main():
    """
    Main function to run the CLI
    """
    try:
        cli()
    except SystemExit as e:
        if e.code != 0:
            raise


if __name__ == "__main__":
    main()
