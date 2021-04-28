import pandas as pd
import pyranges as pr

df_results = pd.read_csv(snakemake.input['results'])

df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)
df_gtf = df_gtf[df_gtf['Feature'] == 'gene']
df_gtf = df_gtf[['Chromosome', 'Start', 'End',
                 'Strand', 'gene_id', 'gene_name', 'gene_biotype']]

genes = set(df_results['geneID'])
df_gtf = df_gtf[df_gtf['gene_id'].isin(genes)]

df_gtf.rename(columns={'gene_id': 'geneID'}) \
      .to_csv(snakemake.output['genes'], index=False)
