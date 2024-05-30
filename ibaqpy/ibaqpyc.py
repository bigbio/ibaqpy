import click

from ibaqpy.commands.features2peptides import features2parquet
from ibaqpy.commands.compute_tpa import tpa_compute
from ibaqpy.commands.datasets_merger import datasets_merge
from ibaqpy.commands.merge_condition_files import merge_condition_generation
from ibaqpy.commands.peptides2proteins import peptides2proteins
from ibaqpy.commands.tsne_visualization import tsne_visualization
import ibaqpy.__init__ as __init__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    version=__init__.__version__, package_name="ibaqpy", message="%(package)s %(version)s"
)

def cli():
    """
    Main entry point for the CLI
    """
    pass

cli.add_command(tpa_compute)
cli.add_command(datasets_merge)
cli.add_command(features2parquet)
cli.add_command(merge_condition_generation)
cli.add_command(peptides2proteins)
cli.add_command(tsne_visualization)


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
