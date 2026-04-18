#!/bin/bash
#SBATCH --job-name=phylo_tree_genus
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=350:00:00
#SBATCH --error ./logs/phylo_tree_genus.e%j
#SBATCH --output ./logs/phylo_tree_genus.o%j


# Define directories
Work_DIR="/mount/britton/Jose/Jobs/Jason/2025-12-15-Phylogram-Tree-Genus-Nextflow-Run/workflow-output/work"  

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
nextflow -c phylogenetic_tree_genus.conf -log ./logs/nextflow.log run phylogenetic_tree_genus.nf -profile MAB -resume --sample Pseudo --type_accession ../data/genus_type_species_type_strain_genus_anchor_assembly_accession.txt

