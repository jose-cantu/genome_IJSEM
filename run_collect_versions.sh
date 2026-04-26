#!/bin/bash 

#SBATCH --job-name=ijsem_versions 
#SBATCH --nodes=1
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output ./logs/collect_versions.%j.out
#SBATCH --error ./logs/collect_versions.%j.err 

cd /mount/britton/Jose/Jobs/Jason/genome_IJSEM
mkdir -p logs docs

export CONDA_EXE="/mount/britton/Jose/conda-envs/mamba/bin/conda"

bash collect_versions_cluster.sh \
  . \
  ../2025-09-22-novel-isolate-closed-genome-trycycler-IJSEM-Criteria-Updates-nextflow-run \
  ../2025-09-30-obtain-ANI-results-type-strain-nextflow-run \
  ../2025-10-04-Obtaining-AAI-Results-and-Tree-Nextflow-run \
  ../2026-04-24-Finalized-Phylogram-Type-Species-Only-Tree-Genus-Nextflow-Run
