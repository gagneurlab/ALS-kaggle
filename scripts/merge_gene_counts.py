import pdb
from pathlib import Path
import pandas as pd


count_dir = Path(snakemake.input['count_dir'])
samples = pd.read_csv(snakemake.input['samples'])['Participant_ID'].tolist()
samples = set(samples).intersection(i.name for i in list(count_dir.iterdir()))
samples = list(samples)

print(f'read {len(samples)} samples')

dfs = list()

for i in samples:
    print(i)
    path = next((count_dir / i).rglob('*.txt'))
    df = pd.read_csv(path, sep='\t', comment='#', index_col=0)
    dfs.append(df[[df.columns[-1]]])

df = pd.concat(dfs, axis=1)

with open(snakemake.output['counts'], 'w') as f:
    f.write('\t'.join(['gene_id'] + samples) + '\n')

df.to_csv(snakemake.output['counts'], mode='a', header=False, sep='\t')
