# Install instructions for VEP v99

The VEP annotation step requires around 2TB of cache files that have to be downloaded beforehand.
Also, it is recommended to run this step on a larger server (> 250GB RAM).

In order to setup VEP with all necessary cache files, please follow the following instructions:

1) install dependencies into some conda environment:
```bash
conda install -c conda-forge glow "pyspark>=3"
conda install -c bioconda "ensembl-vep==99"

auto_install_option="apfc"

echo "n" | vep_install --NO_HTSLIB \
    -v $VERSION -a "$auto_install_option" -g "all" \
    -d ${CONDA_PREFIX}/modules -c "${CACHE_DIR}" \
    -s "homo_sapiens,homo_sapiens_refseq,homo_sapiens_merged" -y "GRCh37" -t \
    -g "CADD,Conservation"

echo "n" | vep_install --NO_HTSLIB \
    -v $VERSION -a "$auto_install_option" -g "all" \
    -d ${CONDA_PREFIX}/modules -c "${CACHE_DIR}" \
    -s "homo_sapiens,homo_sapiens_refseq,homo_sapiens_merged" -y "GRCh38" -t \
    -g "CADD,Conservation"


```

2) download GnomAD:
```bash
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz
```

3) download CADD v1.6:
```bash
#!/bin/bash

for d in v1.6 v1.6/GRCh37 v1.6/GRCh38
do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done


# for GRCh37 / hg19
cd v1.6/GRCh37
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/annotationsGRCh37_v1.6.tar.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi
tar -kzxvf annotationsGRCh37_v1.6.tar.gz

# for GRCh38 / hg38
cd ../../v1.6/GRCh38
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi
tar -kzxvf annotationsGRCh38_v1.6.tar.gz
```

