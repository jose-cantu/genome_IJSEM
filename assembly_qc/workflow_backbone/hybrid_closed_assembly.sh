#!/bin/bash
#SBATCH --job-name=hybrid_closed_assembly
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=120:00:00
#SBATCH --error ./logs/hybrid_closed_assembly.e%j
#SBATCH --output ./logs/hybrid_closed_assembly.o%j


# Define directories
Work_DIR="/mount/britton/Jose/Jobs/Jason/2025-09-22-novel-isolate-closed-genome-trycycler-IJSEM-Criteria-Updates-nextflow-run/workflow-output/work"  

eval "$(conda shell.bash hook)"

# Loading my own conda environment for Nextflow
source activate /mount/britton/Jose/conda-envs/nextflow

export PATH="/mount/britton/Jose/DBs/nextflow_files:$PATH"
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOC=false
export NXF_CONDA_ENABLED=true
export NXF_EXECUTOR=slurm
export NXF_WORK=${Work_DIR}

# Initiate Nextflow job
nextflow -c hybrid_closed_assembly.conf -log ./logs/nextflow.log run hybrid_closed_assembly.nf -profile MAB -resume

