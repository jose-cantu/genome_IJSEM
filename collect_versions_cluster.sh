#!/usr/bin/env bash

OUT="docs/versions.txt"
ENV_OUTDIR="docs/conda_envs"
CONDA_EXE="${CONDA_EXE:-/mount/britton/Jose/conda-envs/mamba/bin/conda}"

mkdir -p docs "$ENV_OUTDIR"
: > "$OUT"

ROOTS=("$@")
if [[ ${#ROOTS[@]} -eq 0 ]]; then
  ROOTS=(".")
fi

section () {
  printf "\n## %s\n" "$1" >> "$OUT"
}

record_env () {
  local label="$1"
  local env_path="$2"
  shift 2

  {
    printf "\n### %s\n" "$label"
    printf "Conda env: %s\n" "$env_path"
    printf "Command: %s run -p %s %s\n" "$CONDA_EXE" "$env_path" "$*"
  } >> "$OUT"

  if [[ -d "$env_path" ]]; then
    "$CONDA_EXE" run -p "$env_path" "$@" >> "$OUT" 2>&1 || true
  else
    echo "MISSING_ENV: $env_path" >> "$OUT"
  fi
}

scan_workflow_files () {
  for root in "${ROOTS[@]}"; do
    [[ -e "$root" ]] || continue
    find "$root" -maxdepth 8 \
      -type d \( -name work -o -name .nextflow \) -prune -o \
      -type f \( -name "*.nf" -o -name "*.conf" -o -name "*.sh" -o -name "*.py" \) \
      -print
  done
}

scan_nextflow_logs () {
  for root in "${ROOTS[@]}"; do
    [[ -e "$root" ]] || continue
    find "$root" -maxdepth 8 \
      -type d \( -name work -o -name .nextflow \) -prune -o \
      -type f \( -name "nextflow.log" -o -name "nextflow.log.*" \) \
      -print
  done
}

section "Run metadata"
{
  echo "Generated: $(date -Is)"
  echo "Host: $(hostname)"
  echo "User: $(whoami)"
  echo "Working directory: $(pwd)"
  echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
  echo "Git branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'not a git repo')"
  echo "Scanned roots:"
  printf '  %s\n' "${ROOTS[@]}"
} >> "$OUT"

section "Nextflow versions from workflow logs"
{
  echo "Nextflow was not executed directly here because the local Nextflow Java runtime is not clean on this node."
  echo "Versions below are parsed from nextflow.log files generated during actual workflow runs."
  echo
} >> "$OUT"

scan_nextflow_logs | sort -u > docs/nextflow_logs_scanned.txt || true

if [[ -s docs/nextflow_logs_scanned.txt ]]; then
  while IFS= read -r log; do
    {
      echo
      echo "### $log"
      grep -m 1 -E "N E X T F L O W.*version|Version:[[:space:]]+[0-9]" "$log" || true
    } >> "$OUT"
  done < docs/nextflow_logs_scanned.txt
else
  echo "No nextflow.log files found in scanned roots." >> "$OUT"
fi

section "Conda environments declared in workflow files"

: > docs/conda_env_paths_used.txt

scan_workflow_files | sort -u > docs/workflow_files_scanned.txt || true

while IFS= read -r wf; do
  grep -hoE "conda[[:space:]]+['\"][^'\"]+['\"]" "$wf" 2>/dev/null \
    | sed -E "s/conda[[:space:]]+['\"]//; s/['\"]$//" \
    >> docs/conda_env_paths_used.txt || true
done < docs/workflow_files_scanned.txt

# Add known envs used in this project, including tools that may be called from shell scripts.
cat >> docs/conda_env_paths_used.txt <<'EOF'
/mount/britton/Jose/conda-envs/nextflow
/mount/britton/Jose/conda-envs/porechop_abi
/mount/britton/Jose/conda-envs/NanoPlot
/mount/britton/Jose/conda-envs/fastp
/mount/britton/Jose/conda-envs/quast
/mount/britton/Jose/conda-envs/checkm
/mount/britton/Jose/conda-envs/checkm2
/mount/britton/Jose/conda-envs/busco
/mount/britton/Jose/conda-envs/yak
/mount/britton/Jose/conda-envs/gtdbtk
/mount/britton/Jose/conda-envs/iqtree
/mount/britton/Jose/conda-envs/ncbi_datasets
/mount/britton/Jose/conda-envs/fastani
/mount/britton/Jose/conda-envs/bakta_env
/mount/britton/Jose/conda-envs/mamba/envs/autocycler_env
/mount/britton/Jose/conda-envs/mamba/envs/dnaapler_env
/mount/britton/Jose/conda-envs/mamba/envs/barrnap_env
/mount/britton/Jose/conda-envs/mamba/envs/ezaai
/mount/britton/Jose/conda-envs/polypolish_polyca_bwa_bwa_mem2_samtools
/mount/britton/Jose/conda-envs/pypolca
EOF

# Remove parameter placeholders and duplicates.
grep -v '\${' docs/conda_env_paths_used.txt | grep -v '^$' | sort -u > docs/conda_env_paths_used.tmp
mv docs/conda_env_paths_used.tmp docs/conda_env_paths_used.txt

while IFS= read -r env_path; do
  [[ -z "$env_path" ]] && continue
  env_name="$(basename "$env_path")"

  {
    echo
    echo "### $env_name"
    echo "Path: $env_path"
  } >> "$OUT"

  if [[ -d "$env_path" ]]; then
    "$CONDA_EXE" list -p "$env_path" > "$ENV_OUTDIR/${env_name}.conda-list.txt" 2>&1 || true
    "$CONDA_EXE" env export -p "$env_path" --no-builds > "$ENV_OUTDIR/${env_name}.environment.yml" 2>&1 || true
    echo "conda list: docs/conda_envs/${env_name}.conda-list.txt" >> "$OUT"
    echo "environment export: docs/conda_envs/${env_name}.environment.yml" >> "$OUT"
  else
    echo "MISSING_ENV" >> "$OUT"
  fi
done < docs/conda_env_paths_used.txt

section "Primary tool versions"

record_env "Porechop_ABI" "/mount/britton/Jose/conda-envs/porechop_abi" bash -lc 'porechop_abi --version || porechop_abi -h | head -20'
record_env "Flye" "/mount/britton/Jose/conda-envs/mamba/envs/autocycler_env" flye --version
record_env "QUAST" "/mount/britton/Jose/conda-envs/quast" quast.py --version
record_env "CheckM" "/mount/britton/Jose/conda-envs/checkm" bash -lc 'checkm version 2>/dev/null || checkm -h | head -20'
record_env "CheckM2" "/mount/britton/Jose/conda-envs/checkm2" checkm2 --version
record_env "BUSCO" "/mount/britton/Jose/conda-envs/busco" busco --version
record_env "Yak" "/mount/britton/Jose/conda-envs/yak" bash -lc 'yak 2>&1 | head -10'
record_env "Barrnap" "/mount/britton/Jose/conda-envs/mamba/envs/barrnap_env" barrnap --version
record_env "NCBI datasets" "/mount/britton/Jose/conda-envs/ncbi_datasets" datasets --version
record_env "GTDB-Tk" "/mount/britton/Jose/conda-envs/gtdbtk" gtdbtk --version
record_env "IQ-TREE 2" "/mount/britton/Jose/conda-envs/iqtree" iqtree2 --version
record_env "FastANI" "/mount/britton/Jose/conda-envs/fastani" fastANI --version
record_env "Bakta" "/mount/britton/Jose/conda-envs/bakta_env" bakta --version
record_env "DNAAPLER" "/mount/britton/Jose/conda-envs/mamba/envs/dnaapler_env" dnaapler --version
record_env "fastp" "/mount/britton/Jose/conda-envs/fastp" fastp --version
record_env "NanoPlot" "/mount/britton/Jose/conda-envs/NanoPlot" NanoPlot --version
record_env "Pypolca" "/mount/britton/Jose/conda-envs/pypolca" pypolca --version
record_env "Polypolish-related env contents" "/mount/britton/Jose/conda-envs/polypolish_polyca_bwa_bwa_mem2_samtools" bash -lc 'which polypolish || true; polypolish --version 2>/dev/null || true; bwa-mem2 version 2>&1 | head -5 || true; samtools --version | head -5 || true'
record_env "EzAAI env contents" "/mount/britton/Jose/conda-envs/mamba/envs/ezaai" bash -lc 'which ezaai || which EzAAI || find "$CONDA_PREFIX" -maxdepth 4 \( -iname "*ezaai*" -o -iname "*.jar" \) | head -40'

section "Databases and web services"
{
  echo
  echo "### GTDB-Tk database"
  echo "GTDBTK_DATA_PATH=/mount/britton/leap/dbs/gtdbtk"

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

section "Workflow files scanned"
cat docs/workflow_files_scanned.txt >> "$OUT"

echo "Wrote $OUT"
echo "Also wrote:"
echo "  docs/conda_env_paths_used.txt"
echo "  docs/nextflow_logs_scanned.txt"
echo "  docs/workflow_files_scanned.txt"
echo "  docs/conda_envs/*.conda-list.txt"
echo "  docs/conda_envs/*.environment.yml"
