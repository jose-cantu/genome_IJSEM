#!/bin/bash
#SBATCH --job-name=ani_scoring_only
#SBATCH --mail-user=jose.cantu@bcm.edu 
#SBATCH --mail-type=END,FAIL 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem=16G 
#SBATCH --time=120:00:00 
#SBATCH --error ./logs/ani_current_workflow.e%j 
#SBATCH --output ./logs/ani_current_workflow.o%j 


# Define directories
Work_DIR="/mount/britton/Jose/Jobs/Jason/2025-09-30-obtain-ANI-results-type-strain-nextflow-run/workflow-output/work"  

eval "$(conda shell.bash hook)"

# Loading my own conda environment for Nextflow
source activate /mount/britton/Jose/conda-envs/nextflow

export NXF_JAVA_HOME="$CONDA_PREFIX"
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false
export NXF_CONDA_ENABLED=true
export NXF_EXECUTOR=local 
export NXF_WORK=${Work_DIR}

# Initiate Nextflow job
nextflow -c ani_current_workflow.conf -log ./logs/nextflow.log run ani_current_workflow.nf -resume -profile standard --sample Pseudo --asm ../data/assembly.polished.fasta --class_summary ../data/gtdbtk.bac120.summary.tsv 

