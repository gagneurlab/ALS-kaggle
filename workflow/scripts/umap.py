import pandas as pd
import matplotlib.pyplot as plt
import umap
import numba
import seaborn as sns
sns.set_theme(style="whitegrid")


numba.set_num_threads(snakemake.threads)

gene_counts = pd.read_table(snakemake.input['counts'], sep='\t')

del gene_counts['gene_id']

X = gene_counts.transpose()
embed_counts = umap.UMAP(random_state=42).fit(X)

df_samples = pd.read_csv(snakemake.input['samples'], sep='\t')
df_samples = df_samples.set_index('Participant_ID').loc[gene_counts.columns]
df_samples['dim1'] = embed_counts.embedding_[:, 0]
df_samples['dim2'] = embed_counts.embedding_[:, 1]

plt.figure(figsize=(5, 5), dpi=250)
annotation_column = snakemake.wildcards.annotation_column
ax = sns.scatterplot(data=df_samples, x="dim1",
                     y="dim2", hue=annotation_column)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(snakemake.output['plot'], bbox_inches='tight')
