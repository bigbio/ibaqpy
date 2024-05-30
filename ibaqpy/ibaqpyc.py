import click

from ibaqpy.commands.compute_tpa import tpa_compute
from ibaqpy.commands.datasets_merger import datasets_merge

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

cli.add_command(tpa_compute)
cli.add_command(datasets_merge)

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