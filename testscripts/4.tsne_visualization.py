# import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


# function to compute principal components
def compute_pca(df, n_components=5) -> pd.DataFrame:
    """
    Compute principal components for a given dataframe.

    Parameters
    ----------
    df : pd.DataFrame
    A dataframe with samples as rows and features as columns.
    n_components : int
    Number of principal components to be computed.

    Returns
    -------
    df_pca : pd.DataFrame
    A dataframe with the principal components.
    """

    pca = PCA(n_components=n_components)
    pca.fit(df)
    df_pca = pca.transform(df)

    df_pca = pd.DataFrame(
        df_pca, index=df.index, columns=[f"PC{i}" for i in range(1, n_components + 1)]
    )

    return df_pca


def compute_tsne(df_pca, n_components=2, perplexity=30, learning_rate=200, n_iter=2000):
    """
    Compute t-SNE components from PCA components.

    This function applies t-SNE (t-Distributed Stochastic Neighbor Embedding) to the input DataFrame,
    which is expected to contain PCA components with samples as rows. The output is another DataFrame
    that contains t-SNE components, also with samples as rows.

    Parameters
    ----------
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

    Returns
    -------
    df_tsne : pandas DataFrame
        Output DataFrame containing t-SNE components. Rows are samples and columns are t-SNE components.

    Example
    -------
     df_pca = pd.DataFrame(data, columns=['PC1', 'PC2', 'PC3'])
     df_tsne = compute_tsne(df_pca)
    """

    tsne = TSNE(n_components=n_components, perplexity=perplexity, learning_rate=learning_rate, n_iter=n_iter)
    tsne_results = tsne.fit_transform(np.asarray(df_pca))

    tsne_cols = [f'tSNE{i + 1}' for i in range(n_components)]

    df_tsne = pd.DataFrame(data=tsne_results, columns=tsne_cols)
    df_tsne.index = df_pca.index
    return df_tsne


def plot_tsne(df, x_col, y_col, hue_col, file_name):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    sns.scatterplot(x=x_col, y=y_col, hue=hue_col, data=df, ax=ax)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f'{x_col} vs {y_col} with {hue_col} information')
    # set legend inside the plot left upper corner
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.subplots_adjust(right=0.8)
    plt.savefig(file_name)


# import data with first column as index
path = "/Users/enrique/projects/local/heart-proteomics"
df = pd.read_table(
    f'{path}/data/processed/heart-Intensities-reshaped-IbaqNorm-KNNImputed-BatchCorrected-OutlierRemoved.tsv',
    index_col=0,
    header=0
)

# generate a t-SNE plot from PCA components
# save the plot as a png file
df_pca = compute_pca(df.T, n_components=15)
df_tsne = compute_tsne(df_pca)

# Add batch info column from sample ids (index)
batch = df_tsne.index.str.split("-").str[0]
df_tsne['batch'] = batch

# plot the t-SNE components tSNE1 vs tSNE2 with batch information using seaborn
plot_tsne(df_tsne, 'tSNE1', 'tSNE2', 'batch', '5.tsne_plot_with_batch_information.pdf')

# import metadata file
df_metadata = pd.read_csv(f'{path}/data/meta_data/heart_samples_meta_data.tsv',
                          sep='\t',
                          index_col=False)
# Set index with column 'sample_id'
df_metadata.set_index('sample_id', inplace=True)

# remove batch column from metadata dataframe
df_metadata.drop('batch', axis=1, inplace=True)

# perform a left join by index between the t-SNE dataframe and the metadata dataframe
df_tsne_metadata = df_tsne.join(df_metadata, how='left')

# keep only entries matching "LV", "RV", "LA", "RA" from the column "organism_part_abbrev"
df_tsne_metadata = df_tsne_metadata[df_tsne_metadata['organism_part_abbrev'].isin(["LV", "RV", "LA", "RA"])]

# plot the t-SNE components tSNE1 vs tSNE2 labeled by metadata using seaborn
plot_tsne(df_tsne_metadata, 'tSNE1', 'tSNE2', 'organism_part_abbrev', '6.tsne_plot_with_metadata_information.pdf')
