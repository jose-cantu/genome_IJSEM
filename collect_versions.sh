#!/usr/bin/env bash

OUT="docs/versions.txt"
mkdir -p docs
: > "$OUT"

section () {
  printf "\n## %s\n" "$1" >> "$OUT"
}

record () {
  local label="$1"
  shift
  {
    printf "\n### %s\n" "$label"
    printf "Command: %s\n" "$*"
    "$@" 2>&1 || true
  } >> "$OUT"
}

record_env () {
  local label="$1"
  local env_path="$2"
  shift 2
  {
    printf "\n### %s\n" "$label"
    printf "Conda env: %s\n" "$env_path"
    printf "Command: conda run -p %s %s\n" "$env_path" "$*"
    conda run -p "$env_path" "$@" 2>&1 || true
  } >> "$OUT"
}

section "Run metadata"
{
  echo "Generated: $(date -Is)"
  echo "Host: $(hostname)"
  echo "User: $(whoami)"
  echo "Working directory: $(pwd)"
  echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
  echo "Git branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'not a git repo')"
} >> "$OUT"

section "Workflow engine"
record "Nextflow" nextflow -version

section "Assembly and QC tools"
record_env "Porechop_ABI" "/mount/britton/Jose/conda-envs/porechop_abi" porechop_abi --version
record_env "Flye" "/mount/britton/Jose/conda-envs/mamba/envs/autocycler_env" flye --version
record_env "QUAST" "/mount/britton/Jose/conda-envs/quast" quast.py --version
record_env "CheckM" "/mount/britton/Jose/conda-envs/checkm" checkm -h
record_env "BUSCO" "/mount/britton/Jose/conda-envs/busco" busco --version
record_env "Yak" "/mount/britton/Jose/conda-envs/yak" yak
record_env "Barrnap" "/mount/britton/Jose/conda-envs/mamba/envs/barrnap_env" barrnap --version

section "Taxonomy and phylogeny tools"
record_env "NCBI datasets" "/mount/britton/Jose/conda-envs/ncbi_datasets" datasets --version
record_env "GTDB-Tk" "/mount/britton/Jose/conda-envs/gtdbtk" gtdbtk --version
record_env "IQ-TREE 2" "/mount/britton/Jose/conda-envs/iqtree" iqtree2 --version

section "AAI tools"
record_env "EzAAI env check" "/mount/britton/Jose/conda-envs/mamba/envs/ezaai" bash -lc 'which ezaai || which EzAAI || find "$CONDA_PREFIX" -maxdepth 3 -iname "*ezaai*" -o -iname "*.jar" | head -20'

section "Web services / visualization"
{
  echo
  echo "### TYGS"
  echo "Service: Type (Strain) Genome Server"
  echo "Result date from report: 2025-09-30"
  echo "Job ID: 1437b825-85b8-4ebb-868c-ab4437c4d4f3"

  echo
  echo "### iTOL"
  echo "Service: Interactive Tree Of Life"
  echo "Used for final tree visualization and manual label editing"
} >> "$OUT"

section "Workflow files inspected"
{
  find . -maxdepth 4 -type f \( -name "*.nf" -o -name "*.conf" -o -name "*.sh" \) | sort
} >> "$OUT"

echo "Wrote $OUT"
