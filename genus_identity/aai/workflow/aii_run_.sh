#!/bin/bash
#SBATCH --job-name=ani_scoring_only
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=120:00:00
#SBATCH --error ./logs/aii_run.e%j
#SBATCH --output ./logs/aii_run.o%j


# Define directories
Work_DIR="/mount/britton/Jose/Jobs/Jason/2025-10-04-Obtaining-AAI-Results-and-Tree-Nextflow-run/workflow-output/work"  

eval "$(conda shell.bash hook)"

# Loading my own conda environment for Nextflow
source activate /mount/britton/Jose/conda-envs/nextflow

export PATH="/mount/britton/Jose/DBs/nextflow_files:$PATH"
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false
export NXF_CONDA_ENABLED=true
export NXF_EXECUTOR=slurm
export NXF_WORK=${Work_DIR}

# Initiate Nextflow job
nextflow -c aii_run_.conf -log ./logs/nextflow.log run aii_run_.nf -profile MAB -resume --sample Pseudo --asm ../data/assembly.polished.fasta --type_accessions ../data/type_accessions.txt 

