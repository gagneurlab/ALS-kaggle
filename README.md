End ALS Kaggle Challenge
==============================

Repository: [GitHub - gagneurlab/ALS](https://github.com/gagneurlab/ALS)

## When the outlier is the signal

Amyotrophic Lateral Sclerosis (ALS) is a neurodegenerative disease for which evidence of molecular causes are not yet well understood. Genome-wide association studies and other case-control approaches have tried to identify genes that could possibly be involved in the development and progression of the disease. However, for most cases, the genetic causes remain a mystery.

We focus on Task 1 of the End ALS kaggle challenge, asking whether there are one or more pathways that could contribute to ALS.

### Approach

We asked whether combining DNA and RNA data could be employed to discover novel genes and biological pathways implicated in ALS. Given the heterogeneity of genetic causes of ALS reported in the literature, we assumed that no single gene or pathway would collectively explain disease cases. Therefore, instead of comparing cases against controls, we searched for potential aberrations that could be specific for a single or a few individuals.

To this end, we first identified expression outliers, which are defined as genes exhibiting aberrant expression levels. We therefore used the machine learning technique of a denoising autoencoder, that we previously adapted to work with RNA-sequencing data [1]. To ensure that the outliers are not technical artefacts of some sort but driven by genetics, we filtered them down to genes that carried rare genetic variants that could be likely to disrupt gene expression. Specifically, we identified and narrowed down the set of variants based on whether they were predicted to have a high or moderate impact or to affect splicing. This analysis yielded a total of 109 genetically-supported outlier genes.

Among these, four genes corresponded to catalogued ALS mutations or known ALS genes, supporting the validity of the approach and possibly yielding a more detailed genetic diagnosis of those patients. 

In order to put potential ALS gene candidates into perspective, we searched for new ALS genes using a two-fold strategy. First, we considered genes that are functionally linked to established ALS genes by applying network diffusion on the STRING gene network. This led to the identification of 16 genetically-supported outliers which are functionally linked to ALS genes or are well established disease genes. Second, we asked whether new potential pathways would emerge as clusters of functionally related genetically-supported outliers. This second analysis revealed gene groups that are members of nucleopores, kinetochores, DNA repair and machinery, as well as ribosomal RNA biogenesis.

Altogether, these results demonstrate the applicability of expression outlier and network diffusion analysis to genetic predispositions of ALS. Our results indicate the involvement of new genes and pathways and supports the view of a complex genetic architecture of ALS. 

*[1] Brechtmann, Felix, Christian Mertes, Agnė Matusevičiūtė, Vicente A. Yépez, Žiga Avsec, Maximilian Herzog, Daniel M. Bader, Holger Prokisch, and Julien Gagneur. 2018. “OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data.” American Journal of Human Genetics 103 (6): 907–17.*

### Results

Our results can be found as a presentation in [`reports/ALS_slides.pdf`](https://github.com/gagneurlab/ALS/tree/master/reports/ALS_slides).

## Reproducibility

All analyses are impemented in a fully reproducible [Snakemake](https://snakemake.github.io) pipeline. In order to execute the pipeline, you need to first install all the depedencies, which you can do using [anaconda](https://docs.conda.io). Next, you need to link the input data to connect it to the workflow.

**Rule graph of the Snakemake workflow**

![rulegraph](reports/figures/dependency/rulegraph.png "Rule graph of the Snakemake workflow")

### Installation

Please install the following conda environment from the repository root.

```commandline
conda env create -f environment.yml
```

This will create a conda environment called `als`.

**Note:** You can alternatively use [mamba](https://mamba.readthedocs.io/en/latest/index.html) for faster install times.

After a fresh install of the environment, you need to install the correct Ensembl release.
This only needs to be done once.

```commandline
conda activate als
pyensembl install --release 84 --species human
```

For VEP, the setup is a bit more complicated.
In case you don't already have a VEP setup, please refer to `references/setup_VEP.md`.

### Prepare the input data

Before executing the pipeline, you need to prepare the input data.
If you prepare the data as described here, you won't need to change any paths in the `configs/config.yml`.

1. Create a link in `data/` called `raw` that links to your location of `end-als` from the [End ALS Kaggle Challenge](https://www.kaggle.com/alsgroup/end-als).

2. Create a directory or link to an empty directory elsewhere in `data/` called `processed`.
   We recommend that use a symlink, as the output can get quite large

3. In `data/external/hg38/` download or copy a annotation GTF from Ensembl and FASTA reference file for the assembly hg38.
   Here, we call the GTF `Homo_sapiens.GRCh38.84.gtf` and the FASTA `GRCh38.primary_assembly.genome.fa`, which you could also modify in the `configs/config.yml`.

4. Please make sure that all necessary inputs are either linked to `data/external/vep` according to `configs/config.yml` or specify different locations of the `CADD`, `LOFTEE` and `gnomAD` directories.

### Run the pipeline

The workflow is implemented as a [Snakemake](https://snakemake.github.io) pipeline.
Snakemake is installed from the conda command above.
Before running the pipeline, please adapt the `configs/config.yaml` file so that it links to the correct input.

After adjusting the paths according to your input locations, you can call the pipeline from the repository root via

```commandline
snakemake -n
```

This should give you an overview of the jobs that Snakemake would want to run.
This is a good point to check of the input and output paths are set correctly.
If this step completes without errors, you can continue to execute the steps shown.
We recommend you to start by visualizing the actual workflow as a directed acyclic graph first, for a better understanding of what is run.
The following command creates a dependency file of the rules (`reports/figures/dependency/rulegraph.png`) as well as all runs (`reports/figures/dependency/dag.png`):

```commandline
snakemake dependency -j 1
```

where `-j 1` allows the pipeline to use at most 1 core.

To run the complete pipeline with e.g. 10 cores, call

```commandline
snakemake -j 10
```

All independent steps will run in parallel if you specify more than 1 core.

Folder Structure
----------------

    ├── LICENSE
    |
    ├── README.md          <- The top-level README for developers using this project.
    |
    ├── configs            <- Config files for Snakemake pipeline
    |
    ├── data
    │   ├── external       <- Data from third party sources. You need to include reference
    |   |                     files in external/hg38/ directory
    |   |
    │   ├── interim        <- Intermediate data that has been transformed.
    |   |
    │   ├── processed      <- The final, canonical data sets for modeling.
    |   |                     This directory will be created autmatically once you run the pipeline.
    |   |                     We recommend that you create a link to an existing directory, as the 
    |   |                     output can be pretty large. 
    |   |
    │   └── raw            <- The original, immutable data dump. You need to create a link
    |                         to the end-als data here
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── workflow           <- Snakemake workflow
    |   ├── notebooks      <- Notebooks to be rendered during pipeline execution
    |   ├── scripts        <- Scripts to be run during pipeline execution
    |   └── Snakefile      <- Snakefile with all the Snakemake rules
    |
    ├── environment.yml   <- Conda environment for reproducing the pipeline
    |
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io
