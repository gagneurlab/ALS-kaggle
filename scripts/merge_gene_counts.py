from pathlib import Path
import pandas as pd


samples = pd.read_csv(snakemake.input['samples'])['Participant_ID'].tolist()
count_dir = Path(snakemake.input['count_dir'])

dfs = list()

for i in samples:
    path = next((count_dir / i).rglob('*.txt'))

    df = pd.read_csv(path, sep='\t', index_col=0)
    dfs.append(df[[df.columns[-1]]])


df = pd.concat(dfs, axis=1)

with open(snakemake.output['counts'], 'w') as f:
    f.write('\t'.join(samples) + '\n')

df.to_csv(snakemake.output['counts'], mode='a', header=False, sep='\t')
