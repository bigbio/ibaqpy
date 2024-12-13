import re

import click
from ibaqpy.ibaq.combiner import Combiner


@click.command("datasets_merge", short_help="Merge ibaq results from compute_ibaq")
@click.option(
    "--data_folder",
    "-d",
    help="Data dolfer contains SDRFs and ibaq CSVs.",
    required=True,
)
@click.option("--output", "-o", help="Output file after batch effect removal.", required=True)
@click.option(
    "--covariate",
    "-c",
    default=None,
    help="Indicator included in covariate consideration when datasets are merged.",
)
@click.option("--organism", help="Organism to keep in input data.", default="HUMAN")
@click.option(
    "--covariate_to_keep",
    "-k",
    help="Keep tissue parts from metadata, e.g. 'LV,RV,LA,RA'.",
    default=None,
)
@click.option(
    "--non_missing_percent_to_keep",
    "-m",
    help="non-missing values in each group.",
    default=0.3,
)
@click.option(
    "--skip_outliers_removal",
    help="Skip removing outliers in all datasets.",
    default=False,
    is_flag=True,
)
@click.option(
    "--n_components",
    help="Number of principal components to be computed.",
    default=None,
)
@click.option("--min_cluster", help="The minimum size of clusters.", default=None)
@click.option(
    "--min_sample_num",
    help="The minimum number of samples in a neighborhood for a point to be considered as a core point.",
    default=None,
)
@click.option("--n_iter", help="Number of iterations to be performed.", default=None)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    help="Output debug information.",
    default=False,
    is_flag=True,
)
@click.pass_context
def datasets_merge(
    ctx,
    data_folder: str,
    output: str,
    covariate: str,
    organism: str,
    covariate_to_keep: list,
    non_missing_percent_to_keep: float,
    skip_outliers_removal: bool,
    n_components: int,
    min_cluster: int,
    min_sample_num: int,
    n_iter: int,
    verbose: bool,
):
    if covariate_to_keep:
        covariate_to_keep = re.split(",\s*", covariate_to_keep)
    combiner = Combiner(data_folder=data_folder, covariate=covariate)
    combiner.imputer(covariate_to_keep)
    if not skip_outliers_removal:
        combiner.outlier_removal(n_components, min_cluster, min_sample_num, n_iter)
    combiner.batch_correction(n_components, covariate_to_keep)
    result = combiner.df_corrected
    result.to_csv(output, sep=",", index=True)
