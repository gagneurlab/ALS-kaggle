import pdb
import pandas as pd
from mmsplice_scripts.data.maf import VariantMafDB
from splicing_outlier_prediction.result import SplicingOutlierResult
from splicing_outlier_prediction.cat_dataloader import CatInference
from count_table import CountTable


results = SplicingOutlierResult.read_csv(snakemake.input['pred'])
results = results[results['tissue'] == snakemake.params['tissue']]

maf = VariantMafDB(snakemake.input['maf'])
results = results.filter_maf(population=maf)

# Filter protein_coding
results.df = results.df[results.df['gene_type'].str.contains('protein_coding')]

# Order gene name
# TO FACTOR: move this to junction_annotation package
results.df['gene_name'] = results.df['gene_name'].map(
    lambda x: ';'.join(sorted(set(x.split(';')))))

results.gene.to_csv(snakemake.output['gene'])
