#!/bin/bash
#SBATCH --job-name=POCPu  
#SBATCH --output=./logs/consprot_%j.out
#SBATCH --error=./logs/consprot_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1 
#SBATCH --mem=64G
#SBATCH --partition=main
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=FAIL,END  # Send email if the job fails or ends

# Load the Anaconda module
module load anaconda3/2023.03-1

# Initialize conda for the shell
source ~/.bashrc

# Activate conda environment for nextflow/consprot 
conda activate /mount/britton/Jose/conda-envs/consprot-dev || { echo "Failed to activate conda environment consprot"; exit 1; }

mkdir -p "$PWD/tmp/consprot_test" 

export TMPDIR="$PWD/tmp/consprot_test" 
unset JAVA_CMD 
export JAVA_HOME="$CONDA_PREFIX/lib/jvm"
export PATH="$CONDA_PREFIX/bin:$JAVA_HOME/bin:$PATH"

export NXF_OPTS="-Djava.io.tmpdir=$TMPDIR -Xms500M -Xmx2G"
export NXF_ANSI_LOG=false
export NXF_CONDA_ENABLED=true

# Diagnostics to make sure everything is running smoothly 
which java 
java -version 
which nextflow 
nextflow -version 
 
# Run consprot workflow directly with the older Nextflow version that compiles this workflow
NF="/share/apps/anaconda3/2023.03-1/bin/nextflow"

echo "Using nextflow for Consprot run:"
echo "$NF"
"$NF" -version 

"$NF" \
    -log "$PWD/logs/nextflow.consprot.log" \
    run /mount/britton/Jose/DBs/github_repos_with_commit_updates_ahead_of_conda_package_tool/consprot/consprot/workflow/main.nf \
    -params-file "$PWD/params.yaml" \
    -c "$PWD/nextflow.config" \
    -resume \
    -w /mount/britton/Jose/Jobs/Jason/2026-05-19-POCPu-Psuedo-Strain-Plus-Panel-From-POCP-Run/workflow-output/work \
    || { echo "Consprot POCPu failed"; exit 1; }

echo "Consprot POCPu classification completed at: $(date)"

