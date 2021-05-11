import pdb
import requests
import pandas as pd


html = requests.get('https://alsod.ac.uk/')
df = pd.read_html(html.text)[1]

df = df.rename(columns={'Gene symbol': 'gene_name', 'Category': 'category'})
del df['Gene name']
del df['Unnamed: 0']
df['cite'] = 'ALSoD'

df.to_csv(snakemake.output['genes'], index=False)
