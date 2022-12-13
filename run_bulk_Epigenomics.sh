#!/bin/bash
cores=$1
params=$2

/bin/bash -c ". /opt/conda/bin/activate"

cd /mnt/
snakemake --snakefile Snakefile_bulk_Epigenomics.py --cores $cores --forcerun all ${params}
