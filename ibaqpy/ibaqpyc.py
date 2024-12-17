import click
from ibaqpy.commands.features2peptides import features2parquet
from ibaqpy.commands.peptides2protein import peptides2protein
from ibaqpy.commands.tsne_visualization import tsne_visualization
import ibaqpy.__init__ as __init__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    version=__init__.__version__,
    package_name="ibaqpy",
    message="%(package)s %(version)s",
)
def cli():
    """
    Main entry point for the CLI
    """

cli.add_command(features2parquet)
cli.add_command(peptides2protein)
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
