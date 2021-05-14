# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian3]
#     language: python
#     name: conda-env-anaconda-florian3-py
# ---

# +
import os
import pandas as pd

import json
import yaml

import pyspark
from pyspark.sql import SparkSession
import pyspark.sql.types as t
import pyspark.sql.functions as f

import glow
# -


snakefile_path = os.getcwd() + "/../Snakefile"
snakefile_path

try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args

    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'annotate_variants_with_vep',
        default_wildcards={
        }
    )

print(json.dumps(snakemake.__dict__, indent=2))

# +
MEM = os.popen("ulimit -m").read()
if MEM.startswith("unlimited"):
    print("Memory not constrained, using all available memory...")
    import psutil
    MEM = psutil.virtual_memory().available / 1024
MEM = int(MEM)

N_CPU = int(os.popen("nproc").read())

print("memory: %dk" % MEM)
print("number of cores: %d" % N_CPU)

# +
# MEM = int(MEM * 0.8)

# +
os.environ['PYSPARK_SUBMIT_ARGS'] = " ".join([
    '--driver-memory %dk' % MEM,
    'pyspark-shell'
])
os.environ['PYSPARK_SUBMIT_ARGS']

MAX_FAILURES=4

spark = (
    SparkSession.builder
    .appName('desmi_inject_gnomad')
    .config("spark.jars.packages", ",".join([
        "io.projectglow:glow-spark3_2.12:1.0.0",
    ]))
    .config("spark.local.dir", os.environ.get("TMP"))
    .config("spark.master", f"local[{N_CPU},{MAX_FAILURES}]")
    .config("spark.sql.shuffle.partitions", "2001")
    .config("spark.sql.execution.arrow.enabled", "true")
    .config("spark.driver.maxResultSize", "48G")
    .config("spark.task.maxFailures", MAX_FAILURES)
    .getOrCreate()
)
glow.register(spark)
spark


# +
# INPUT_VCF  = '/s/raw/als/kaggle/end-als/genomics-data/AnswerALS_subset_annovar.hg38_anno_and_geno.no_intergenic.vcf.gz'
# OUTPUT_PATH = '/s/project/kaggle-als/vep_annotations/'
# OUTPUT_PQ = OUTPUT_PATH + "/ensembl_bychr.parquet"
# -

INPUT_VCF  = snakemake.input["vcf"]
INPUT_VCF

OUTPUT_PQ = snakemake.output["vep"]
OUTPUT_PQ

# +
# HUMAN_GENOME_VERSION="hg19"
# ASSEMBLY="GRCh37"
HUMAN_GENOME_VERSION="hg38"
ASSEMBLY="GRCh38"

GNOMAD_VCF=f"/s/raw/gnomad/2.1.1/{HUMAN_GENOME_VERSION}/gnomad.genomes.r2.1.1.sites.vcf.gz"
FASTA=snakemake.input["fasta"]
GTF=snakemake.input["gtf"]
VEP_CMD = snakemake.params["vep_annotate_cmd"]

VEP_CACHE="/opt/modules/i12g/conda-ensembl-vep/99/cachedir"
VEP_PLUGIN_DIR=f"{VEP_CACHE}/Plugins/99_GRCh38"

LOFTEE_DATA_DIR=f"/s/raw/loftee/{ASSEMBLY}"
LOFTEE_PATH=f"{VEP_PLUGIN_DIR}/loftee"
MAXENTSCAN_DATA_DIR=f"{LOFTEE_PATH}/maxEntScan"

CADD_DIR=f"/s/raw/cadd/v1.6/{ASSEMBLY}"
CADD_WGS_SNV=f"{CADD_DIR}/whole_genome_SNVs.tsv.gz"
CADD_INDEL={
    "GRCh37": f"{CADD_DIR}/InDels.tsv.gz",
    "GRCh38": f"{CADD_DIR}/gnomad.genomes.r3.0.indel.tsv.gz",
}[ASSEMBLY]

vep_cmd=" ".join([
#     "/opt/modules/i12g/anaconda/3-2019.10/bin/conda run --prefix /opt/modules/i12g/conda-ensembl-vep/99",
#     "vep",
    VEP_CMD,
    "--output_file STDOUT",
    "--format vcf",
    f"--cache --offline --dir={VEP_CACHE}",
    "--force_overwrite",
    "--no_stats",
    "--tab",
    "--merged",
    f"--assembly {ASSEMBLY}",
    f"--fasta {FASTA}",
    "--species homo_sapiens",
    "--everything",
    "--allele_number",
    "--total_length",
    "--numbers",
    "--symbol",
    "--hgvs",
    "--ccds",
    "--uniprot",
    "--af",
    "--af_1kg",
    "--af_esp",
    "--af_gnomad",
    "--max_af",
    "--pubmed",
    "--canonical",
    "--biotype",
    "--sift b",
    "--polyphen b",
    "--appris",
    "--domains",
    "--protein",
    "--regulatory",
    "--tsl",
# Does not really work with hg38
#     "--plugin " + ",".join([
#         "LoF",
#         f"loftee_path:{LOFTEE_PATH}",
#         f"gerp_bigwig:{LOFTEE_DATA_DIR}/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
#         f"human_ancestor_fa:{LOFTEE_DATA_DIR}/human_ancestor.fa.gz",
#         f"conservation_file:{LOFTEE_DATA_DIR}/loftee.sql",
#     ]),
#     f"--plugin LoF,human_ancestor_fa:{LOFTEE_DATA_DIR}/human_ancestor.fa.gz,conservation_file:{LOFTEE_DATA_DIR}/loftee.sql,loftee_path:{LOFTEE_PATH}",
    "--plugin Condel",
    f"--plugin MaxEntScan,{MAXENTSCAN_DATA_DIR}",
    "--plugin Blosum62",
    "--plugin miRNA",
    f"--plugin CADD,{CADD_WGS_SNV},{CADD_INDEL}",
])
#    --gtf $GTF
#    --plugin LoFtool
#    --plugin dbscSNV
#    --custom ${GNOMAD_VCF},Gnomad_2.1.1,vcf,overlap
#    --plugin AncestralAllele

print(vep_cmd)
# -

df = (
    spark
    .read
    .option("flattenInfoFields", False)
    .format('vcf')
    .load(INPUT_VCF)
)

df = glow.transform("split_multiallelics", df)

df.printSchema()

df = df.withColumn("names", f.array([f.concat(
    f.col('contigName'),
    f.lit(":"),
    f.col('start') + 1,
    f.lit(":"),
    f.col('referenceAllele'),
    f.lit(">"),
    f.col('alternateAlleles')[0]
)]))

df.limit(10).toPandas()

# +
import json
import shlex

input_df = df.select([
    f.col('contigName'),
    f.col('start'),
    f.col('end'),
    f.col('names'),
    f.col('referenceAllele'),
    f.col('alternateAlleles'),
])

vep_transformed_df = glow.transform(
    "pipe",
    input_df,
#     cmd=json.dumps(shlex.split("cat | grep -v '^##'")),
    cmd=json.dumps(shlex.split(vep_cmd)),
    inputFormatter='vcf',
    inVcfHeader='infer',
    outputFormatter='csv',
#     outQuote="##",
    outHeader=True,
    outDelimiter="\t",
)
# -

vep_transformed_df.printSchema()

parsed_df = (
    vep_transformed_df
    .withColumn("#Uploaded_variation", f.split(f.col("#Uploaded_variation"), '[:>]'))
    .select([
        f.col("#Uploaded_variation")[0].alias("chrom"),
        f.col("#Uploaded_variation")[1].alias("pos"),
        f.col("#Uploaded_variation")[2].alias("ref"),
        f.col("#Uploaded_variation")[3].alias("alt"),
        "*",
    ])
    .withColumn("ref", f.when(f.col("ref") == '-', "").otherwise(f.col("ref")))
    .withColumn("alt", f.when(f.col("alt") == '-', "").otherwise(f.col("alt")))
    .withColumn("CANONICAL", f.when(f.col("CANONICAL") == 'YES', True).otherwise(False))
    .withColumn("condel_prediction", f.regexp_replace(f.col("Condel"), "\\(.*", ""))
    .withColumn("condel_score", f.regexp_replace(f.col("Condel"), "(.*\\()|\\)", ""))
    .withColumn("sift_prediction", f.regexp_replace(f.col("SIFT"), "\\(.*", ""))
    .withColumn("sift_score", f.regexp_replace(f.col("SIFT"), "(.*\\()|\\)", ""))
    .withColumn("polyphen_prediction", f.regexp_replace(f.col("PolyPhen"), "\\(.*", ""))
    .withColumn("polyphen_score", f.regexp_replace(f.col("PolyPhen"), "(.*\\()|\\)", ""))
    .drop(
        "#Uploaded_variation", 
        "Location", 
        "Allele",
        "Condel",
        "SIFT",
        "PolyPhen",
    )
)
parsed_df.printSchema()

# +
# x = parsed_df.filter(f.col("Condel") != '-').limit(10).toPandas()
# x
# -

dtypes = {
    "consequence": t.ArrayType(t.StringType()),
    "Existing_variation": t.ArrayType(t.StringType()),
    "ALLELE_NUM": t.IntegerType(),
    "DISTANCE": t.IntegerType(),
    "STRAND": t.ShortType(),
    "FLAGS": t.ArrayType(t.StringType()),
    "HGNC_ID": t.IntegerType(),
#     "CANONICAL": t.BooleanType(), # needs manual check if column equals "CANONICAL"
    "TREMBL": t.ArrayType(t.StringType()),
    "REFSEQ_MATCH": t.ArrayType(t.StringType()),
    "GENE_PHENO": t.BooleanType(),
    "sift_score": t.FloatType(),
    "polyphen_score": t.FloatType(),
    "EXON": t.ArrayType(t.IntegerType()),
    "INTRON": t.ArrayType(t.IntegerType()),
    "HGVS_OFFSET": t.IntegerType(),
    "AF": t.ArrayType(t.FloatType()),
    "AFR_AF": t.ArrayType(t.FloatType()),
    "AMR_AF": t.ArrayType(t.FloatType()),
    "EAS_AF": t.ArrayType(t.FloatType()),
    "EUR_AF": t.ArrayType(t.FloatType()),
    "SAS_AF": t.ArrayType(t.FloatType()),
    "AA_AF": t.ArrayType(t.FloatType()),
    "EA_AF": t.ArrayType(t.FloatType()),
    "gnomAD_AF": t.ArrayType(t.FloatType()),
    "gnomAD_AFR_AF": t.ArrayType(t.FloatType()),
    "gnomAD_AMR_AF": t.ArrayType(t.FloatType()),
    "gnomAD_ASJ_AF": t.ArrayType(t.FloatType()),
    "gnomAD_EAS_AF": t.ArrayType(t.FloatType()),
    "gnomAD_FIN_AF": t.ArrayType(t.FloatType()),
    "gnomAD_NFE_AF": t.ArrayType(t.FloatType()),
    "gnomAD_OTH_AF": t.ArrayType(t.FloatType()),
    "gnomAD_SAS_AF": t.ArrayType(t.FloatType()),
    "MAX_AF": t.FloatType(),
    "MAX_AF_POPS": t.ArrayType(t.StringType()),
    "PUBMED": t.ArrayType(t.StringType()),
    "MOTIF_POS": t.IntegerType(),
    "MOTIF_SCORE_CHANGE": t.FloatType(),
    "Condel": t.StringType(),
    "condel_score": t.FloatType(),
    "condel_prediction": t.StringType(),
    "BLOSUM62": t.IntegerType(),
    "LoF": t.StringType(),
    "LoF_filter": t.StringType(),
    "LoF_flags": t.StringType(),
    "MaxEntScan_ref": t.FloatType(),
    "MaxEntScan_alt": t.FloatType(),
    "MaxEntScan_diff": t.FloatType(),
    "CADD_PHRED": t.FloatType(),
    "CADD_RAW": t.FloatType(),
}

needsMinVal = {
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "AA_AF",
    "EA_AF",
    "gnomAD_AF",
    "gnomAD_AFR_AF",
    "gnomAD_AMR_AF",
    "gnomAD_ASJ_AF",
    "gnomAD_EAS_AF",
    "gnomAD_FIN_AF",
    "gnomAD_NFE_AF",
    "gnomAD_OTH_AF",
    "gnomAD_SAS_AF",
}


def parse_col(name):
    col = f.col(name)
    col = f.when(col == '-', None).otherwise(col)
    if name in dtypes:
        dtype = dtypes[name]
        if isinstance(dtype, t.ArrayType):
            col = f.split(col, ",")
        col = col.cast(dtypes[name])
    if name in needsMinVal:
        col = f.array_min(col)
    
    return col.alias(name)


parsed_df = parsed_df.select(*[
    parse_col(x) for x in parsed_df.columns
])

parsed_df.printSchema()

parsed_df.filter(~ f.isnull("condel_score")).limit(10).toPandas()

# +
# (
#     parsed_df
#     .filter(f.col("Gene").startswith("ENSG"))
#     .sortWithinPartitions(["chrom", "pos"])
#     .write
#     .mode("overwrite")
#     .partitionBy("Gene")
#     .parquet(OUTPUT_PATH + "/ensembl_genes.parquet")
# )
# -

(
    parsed_df
    .filter(f.col("Gene").startswith("ENSG"))
    .sortWithinPartitions(["chrom", "Gene", "pos"]) # hierarchical sorting
    .write
    .mode("overwrite")
    .partitionBy("chrom")
    .parquet(OUTPUT_PQ)
)


