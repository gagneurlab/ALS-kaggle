#!/bin/bash

# source /opt/modules/i12g/anaconda/3-2019.10/bin/activate /opt/modules/i12g/conda-ensembl-vep/99

vep "$@" | grep -v '^##'
