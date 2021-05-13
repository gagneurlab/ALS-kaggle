als
==============================

A short description of the project.

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>


Installation
------------

Please install the following conda environment from the repository root.

```console
conda env create -f environment.yml
```

**Note:** You can alternatively use [mamba](https://mamba.readthedocs.io/en/latest/index.html) for faster install times, if you prefer.

[TODO: notebooks working in pipeline]


Run the pipeline
----------------

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