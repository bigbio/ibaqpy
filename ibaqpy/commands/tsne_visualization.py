import glob
import logging

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from ibaqpy.ibaq.ibaqpy_commons import PROTEIN_NAME, SAMPLE_ID, IBAQ_LOG


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# function to compute principal components
def compute_pca(df, n_components=5) -> pd.DataFrame:
    """
    Compute principal components for a given dataframe.

    Parameters:
    df : pd.DataFrame
        Input dataframe with samples as rows and features as columns.
    n_components : int, optional
        The number of principal components to compute (default is 5).
    Returns:
        df_pca : pd.DataFrame
        A dataframe with the principal components.
    """

    pca = PCA(n_components=n_components)
    pca.fit(df)
    df_pca = pca.transform(df)

    df_pca = pd.DataFrame(
        df_pca, index=df.index, columns=[f"PC{i}" for i in range(1, n_components + 1)]
    )

    plt.rcParams["figure.figsize"] = (12, 6)

    fig, ax = plt.subplots()
    xi = np.arange(1, n_components + 1, step=1)
    y = np.cumsum(pca.explained_variance_ratio_)

    plt.ylim(0.0, 1.1)
    plt.plot(xi, y, marker="o", linestyle="--", color="b")

    plt.xlabel("Number of Components")
    plt.xticks(
        np.arange(0, n_components, step=1)
    )  # change from 0-based array index to 1-based human-readable label
    plt.ylabel("Cumulative variance (%)")
    plt.title("The number of components needed to explain variance")

    plt.axhline(y=0.95, color="r", linestyle="-")
    plt.text(0.5, 0.85, "95% cut-off threshold", color="red", fontsize=16)

    ax.grid(axis="x")
    plt.show()

    return df_pca


def compute_tsne(df_pca, n_components=2, perplexity=30, learning_rate=200, n_iter=2000):
    """
    Compute t-SNE components from PCA components.

    This function applies t-SNE (t-Distributed Stochastic Neighbor Embedding) to the input DataFrame,
    which is expected to contain PCA components with samples as rows. The output is another DataFrame
    that contains t-SNE components, also with samples as rows.

    Parameters:
    df_pca : pandas DataFrame
        Input DataFrame containing PCA components. Rows are samples and columns are PCA components.
    n_components : int, optional
        The number of dimensions for the t-SNE components (default is 2).
    perplexity : float, optional
        The perplexity parameter for t-SNE, which can influence the balance between maintaining
        the local and global structure of the data (default is 30).
    learning_rate : float, optional
        The learning rate for t-SNE (default is 200).
    n_iter : int, optional
        The number of iterations for t-SNE optimization (default is 2000).

    Returns:
    df_tsne : pandas DataFrame
        Output DataFrame containing t-SNE components. Rows are samples and columns are t-SNE components.

    Example
    -------
     df_pca = pd.DataFrame(data, columns=['PC1', 'PC2', 'PC3'])
     df_tsne = compute_tsne(df_pca)
    """

    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        learning_rate=learning_rate,
        n_iter=n_iter,
    )
    tsne_results = tsne.fit_transform(np.asarray(df_pca))

    tsne_cols = [f"tSNE{i + 1}" for i in range(n_components)]

    df_tsne = pd.DataFrame(data=tsne_results, columns=tsne_cols)
    df_tsne.index = df_pca.index
    return df_tsne


def plot_tsne(df, x_col, y_col, hue_col, file_name):
    """
    Generate and save a t-SNE scatter plot from a DataFrame.

    This function creates a scatter plot using seaborn's scatterplot function,
    with the specified columns for the x-axis, y-axis, and hue. The plot is
    customized with labels, a title, and a legend positioned inside the plot.
    The resulting plot is saved to the specified file.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data to plot.
        x_col (str): The column name for the x-axis values.
        y_col (str): The column name for the y-axis values.
        hue_col (str): The column name for the hue (color) values.
        file_name (str): The file path where the plot image will be saved.
    """
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    sns.scatterplot(x=x_col, y=y_col, hue=hue_col, data=df, ax=ax, markers=["o", "+", "x"])
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{x_col} vs {y_col} with {hue_col} information")
    # set legend inside the plot left an upper corner
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=8)
    plt.subplots_adjust(right=0.8)
    plt.savefig(file_name)


@click.command()
@click.option("-f", "--folder", help="Folder that contains all the protein files", required=True)
@click.option(
    "-o",
    "--pattern",
    help="Protein file pattern",
    # TODO: I think we should use instead of pattern the structure of quantms.io for absolute quantification
    required=False,
    default="proteins.tsv",
)
def tsne_visualization(folder: str, pattern: str):
    """
    Generate a t-SNE visualization for protein data from specified files.

    This command-line tool reads protein data files from a specified folder,
    applies PCA and t-SNE for dimensionality reduction, and generates a scatter
    plot of the t-SNE components. The plot is saved as a PDF file.

    Parameters:
        folder (str): The folder containing protein data files.
        pattern (str): The file pattern to match protein files. Defaults to 'proteins.tsv'.
    """
    # get all the files in the folder
    files = glob.glob(f"{folder}/*{pattern}")

    # get the files into pandas selected columns
    # (Proteins accession, Sample ID, Reanalysis accession, Intensity)

    dfs = []  # list of dataframes

    for f in files:
        reanalysis = (f.split("/")[-1].split("_")[0]).replace("-proteins.tsv", "")
        dfs += [
            pd.read_csv(f, usecols=[PROTEIN_NAME, SAMPLE_ID, IBAQ_LOG], sep=",").assign(
                reanalysis=reanalysis
            )
        ]

    total_proteins = pd.concat(dfs, ignore_index=True)

    normalize_df = pd.pivot_table(
        total_proteins,
        index=[SAMPLE_ID, "reanalysis"],
        columns=PROTEIN_NAME,
        values=IBAQ_LOG,
    )
    normalize_df = normalize_df.fillna(0)
    df_pca = compute_pca(normalize_df, n_components=30)
    df_tsne = compute_tsne(df_pca)

    batch = df_tsne.index.get_level_values("reanalysis").tolist()
    df_tsne["batch"] = batch

    # plot the t-SNE components tSNE1 vs tSNE2 with batch information using seaborn
    plot_tsne(df_tsne, "tSNE1", "tSNE2", "batch", "5.tsne_plot_with_batch_information.pdf")

    logger.info(total_proteins.shape)


if __name__ == "__main__":
    tsne_visualization()
