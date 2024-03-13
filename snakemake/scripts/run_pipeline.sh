#!/bin/bash

# Define the path to the Conda environments as a constant
CONDA_ENV_PATH="/lustre/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs"

# Generate a timestamp
timestamp=$(date +"%Y%m%d%H%M")

# Run Snakemake and redirect stdout and stderr to separate log files with the timestamp
snakemake --cores 12 --use-conda --rerun-incomplete --resources gpu=1 --conda-prefix $CONDA_ENV_PATH > "../log/output_${timestamp}_stdout.log" 2> "../log/output_${timestamp}_stderr.log"

