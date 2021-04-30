from mmsplice import predict_save, MMSplice
from mmsplice.vcf_dataloader import SplicingVCFDataloader

dl = SplicingVCFDataloader('grch38', snakemake.input['fasta'],
                           snakemake.input['vcf'])
model = MMSplice()
predict_save(model, dl, output_csv=snakemake.output['result'])

