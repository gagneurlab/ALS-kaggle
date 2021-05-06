import pandas as pd
import matplotlib.pyplot as plt
import umap
import umap.plot
import numba

numba.set_num_threads(snakemake.threads)

annotation_column = snakemake.wildcards.annotation_column
gene_counts = pd.read_table(snakemake.input['counts'], sep = '\t')
del gene_counts['gene_id']

sample_annotation = pd.read_table(
    snakemake.input['sample_annotation'],
    sep = '\t',
    index_col='Participant_ID'
)
labels = sample_annotation.loc[gene_counts.columns][annotation_column].astype(str)

X = gene_counts.transpose()
embed_counts = umap.UMAP().fit(X)

umap.plot.points(
    embed_counts,
    labels=labels,
    theme='fire'
)
plt.savefig(snakemake.output['plot'])
