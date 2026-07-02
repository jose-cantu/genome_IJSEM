#!/bin/bash
#SBATCH --job-name=multipass
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=120:00:00
#SBATCH --error ./logs/multipass.e%j
#SBATCH --output ./logs/multipass.o%j


# Define directories
Work_DIR="/mount/britton/Jose/Jobs/Jason/2026-05-19-BAKTA-Annotation-Pseudo-Strain/workflow-output/work"    #~1

# Loading my own conda environment for Nextflow
source activate /mount/britton/Jose/conda-envs/nextflow

export PATH="/mount/britton/Jose/DBs/nextflow_files:$PATH"
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOC=false
export NXF_CONDA_ENABLED=true
export NXF_EXECUTOR=slurm
export NXF_WORK=${Work_DIR}

# Initiate Nextflow job
nextflow -c multipass.conf -log ./logs/nextflow.log run multipass.nf -profile MAB -resume

