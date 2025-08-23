#!/bin/bash
# This script runs the Snakemake pipeline with the specified configuration

cd "/path/to/your/project/"
snakemake --snakefile "/path/to/your/project/nanoscope/workflow/Snakefile_preprocess.smk" \
          --cores 16 \
          --jobs 100 \
          -p \
          --use-conda \
          --configfile "/path/to/your/project/your_config.yaml"
          