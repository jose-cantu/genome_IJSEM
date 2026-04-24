#!/bin/bash
#SBATCH --job-name=closable
#SBATCH --output=/mount/britton/Jose/Jobs/Jason/2025-12-15-In-silico-filtering-last-effort-in-closing-genome-bash-run/workflow/logs/closable_%j.out 
#SBATCH --error=/mount/britton/Jose/Jobs/Jason/2025-12-15-In-silico-filtering-last-effort-in-closing-genome-bash-run/workflow/logs/closable_%j.err 
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --partition=main
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=FAIL,END  # Send email if the job fails or ends

# Load the Anaconda module
module load anaconda3/2023.03-1

# Initialize conda for the shell
eval "$(conda shell.bash hook)" 

set -euo pipefail                                                    # abort on error/undef var/pipe failure   # control-flow safety

# Defining the output directory for easier reference downstream 
OUT=../workflow-output
# Definig other variables 
T=${SLURM_CPUS_PER_TASK:-8}

ASM="$OUT/pilon/assembly.clean.fasta"

ONT="$OUT/reads/Pseudo.ont.trim.fastq.gz"

mkdir -p "$OUT/closable"/{nanoplot,map,bridges}

ts(){ date -Iseconds; } 

echo -e "\nstarting: NanoPlot [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/NanoPlot                           
if [ ! -s "$OUT/closable/nanoplot/NanoPlot-report.html" ]; then                  
  NanoPlot --fastq "$ONT" -o "$OUT/closable/nanoplot" -t "$T"                     
fi
conda deactivate
echo -e "ending:   NanoPlot [$(ts)]\n---------------------------"

echo -e "\nstarting: ONT coverage quick check [$(ts)]\n---------------------------------"
if [ ! -s "$OUT/closable/coverage.txt" ]; then
  zcat "$ONT" | awk -v gs=2575873 'NR%4==2{b+=length($0)} END{printf "ONT_bases=%d\tcov=%.2fx\n", b, b/gs}' \
    > "$OUT/closable/coverage.txt"                                               # cov ≈ total_bases / genome_size
fi
echo -e "ending:   coverage quick check [$(ts)]\n---------------------------"

echo -e "\nstarting: minimap2 map-ont [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/mamba/envs/polypolish_polyca_bwa_bwa_mem2_bowtie2_samtools_pilon_minimap2
if [ ! -s "$OUT/closable/map/ont2asm.bam.bai" ]; then
  # making sure minimap2 does not emit secondary alignments 
  minimap2 -ax map-ont --secondary=no -t "$T" "$ASM" "$ONT" \
    | samtools sort -@ "$T" -o "$OUT/closable/map/ont2asm.bam" -                  # BAM from stream
  samtools index "$OUT/closable/map/ont2asm.bam"                                  # .bai
fi
echo -e "ending:   minimap2 map-ont [$(ts)]\n---------------------------"
[ -s "$OUT/closable/map/ont2asm.bam.bai" ] || { echo "BAM build failed"; exit 1; }  # guard invariant  # prevents empty read lists

# Compute contig sizes and run one bridging test (top two contigs)

# Compute contig sizes and run multi-window, dual-orientation bridge tests (top two contigs)
echo -e "\nstarting: bridge test (top 2 contigs: multi-window/orientation) [$(ts)]\n---------------------------------"

# Contig name->length table (token = header up to first space)
awk '
  /^>/{ if(NR>1) printf "%s\t%d\n", name, len; name=substr($0,2); sub(/ .*/,"",name); len=0; next }
  {len+=length($0)}
  END{ if(name!="") printf "%s\t%d\n", name, len}
' "$ASM" | sort -k2,2nr > "$OUT/closable/bridges/contigs.tsv"

CTG1=$(awk 'NR==1{print $1}' "$OUT/closable/bridges/contigs.tsv")   # largest contig
LEN1=$(awk 'NR==1{print $2}' "$OUT/closable/bridges/contigs.tsv")
CTG2=$(awk 'NR==2{print $1}' "$OUT/closable/bridges/contigs.tsv")   # second largest
LEN2=$(awk 'NR==2{print $2}' "$OUT/closable/bridges/contigs.tsv")

BAM="$OUT/closable/map/ont2asm.bam"
[ -s "$BAM.bai" ] || { echo "ERROR: $BAM.bai missing -> mapping failed"; exit 1; }   # assert artifact

# Header (add COMBO column)
printf "CTG1\tLEN1\tCTG2\tLEN2\tWIN\tCOMBO\tSHARED_READS\n" > "$OUT/closable/bridges/summary.tsv"

# COMBO legend: TH=CTG1_tail->CTG2_head, TT=tail->tail, HH=head->head, HT=head->tail
for WIN in 20000 50000 100000; do
  # define head/tail intervals for this WIN (cap to contig length)
  TAIL1_START=$(( LEN1>WIN ? LEN1-WIN+1 : 1 ))
  HEAD1_END=$(( LEN1>WIN ? WIN : LEN1 ))
  TAIL2_START=$(( LEN2>WIN ? LEN2-WIN+1 : 1 ))
  HEAD2_END=$(( LEN2>WIN ? WIN : LEN2 ))

  for COMBO in TH TT HH HT; do
    case "$COMBO" in
      TH) R1="${CTG1}:${TAIL1_START}-${LEN1}"; R2="${CTG2}:1-${HEAD2_END}";;
      TT) R1="${CTG1}:${TAIL1_START}-${LEN1}"; R2="${CTG2}:${TAIL2_START}-${LEN2}";;
      HH) R1="${CTG1}:1-${HEAD1_END}";        R2="${CTG2}:1-${HEAD2_END}";;
      HT) R1="${CTG1}:1-${HEAD1_END}";        R2="${CTG2}:${TAIL2_START}-${LEN2}";;
    esac

    # Count shared read IDs (no temp files) filtering by MAPQ 20 when counting shared read IDs 
    BRIDGES=$(comm -12 \
      <(samtools view -q 20 "$BAM" "$R1" | awk '{print $1}' | sort -u) \
      <(samtools view -q 20 "$BAM" "$R2" | awk '{print $1}' | sort -u) | wc -l)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$CTG1" "$LEN1" "$CTG2" "$LEN2" "$WIN" "$COMBO" "$BRIDGES" \
      | tee -a "$OUT/closable/bridges/summary.tsv"
  done
done

echo -e "ending:   bridge test (top 2 contigs) [$(ts)]\n---------------------------"
conda deactivate

echo -e "\nstarting raven [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/mamba/envs/autocycler_env   
if command -v raven >/dev/null 2>&1; then
  mkdir -p "$OUT/raven"
  if [ ! -s "$OUT/raven/raven.fasta" ]; then
    raven --threads "$T" "$ONT" > "$OUT/raven/raven.fasta"                # quick alt assembly
  fi
else
  echo "WARN: raven not found in PATH; skipping" >&2
fi
conda deactivate
echo -e "ending:   raven [$(ts)]\n---------------------------"

echo -e "\nstarting: quast_raven [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/quast
if [ -s "$OUT/raven/raven.fasta" ]; then
  [ -s "$OUT/quast_raven/report.txt" ] || quast -o "$OUT/quast_raven" -t "$T" "$OUT/raven/raven.fasta"
else
  echo "SKIP: $OUT/raven/raven.fasta missing; quast_raven skipped" >&2
fi
conda deactivate
echo -e "ending:   quast_raven [$(ts)]\n---------------------------"


echo -e "\nstarting: ragtag patch  [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/mamba/envs/RagTag

# sanity: inputs readable
[ -r "$OUT/flye/assembly.fasta" ] || { echo "ERR: $OUT/flye/assembly.fasta not readable"; ls -l "$OUT/flye/assembly.fasta"; exit 1; }
[ -r "$OUT/raven/raven.fasta" ] || { echo "ERR: $OUT/raven/raven.fasta not readable"; ls -l "$OUT/raven/raven.fasta"; exit 1; }

# idempotent: skip if patched FASTA already exists
if [ ! -s "$OUT/ragtag_patch/ragtag.patch.fasta" ]; then
  rm -rf "$OUT/ragtag_patch.tmp"                                                          # clean temp dir (safe overwrite)
  mkdir -p "$OUT/ragtag_patch.tmp"
  ragtag.py patch -o "$OUT/ragtag_patch.tmp" -t "$T" --aligner minimap2 --mm2-params "-x asm5" \
    "$OUT/flye/assembly.fasta" "$OUT/raven/raven.fasta"                                   # patch Flye using Raven
  rm -rf "$OUT/ragtag_patch"
  mv "$OUT/ragtag_patch.tmp" "$OUT/ragtag_patch"                                          # promote temp -> final dir
fi

conda deactivate
echo -e "ending:   ragtag patch  [$(ts)]\n---------------------------"



# Validate
echo -e "\nstarting: validate patched [$(ts)]\n---------------------------------"

# 6a) QUAST in its own env
conda activate /mount/britton/Jose/conda-envs/quast
if [ -s "$OUT/ragtag_patch/ragtag.patch.fasta" ]; then
  FINAL_ASM="$OUT/ragtag_patch/ragtag.patch.fasta"                            # set candidate
  [ -s "$OUT/quast_patched/report.txt" ] || quast -o "$OUT/quast_patched" -t "$T" "$FINAL_ASM"
else
  echo "SKIP: patched assembly missing; using $ASM as FINAL_ASM" >&2
  FINAL_ASM="$ASM"
fi
conda deactivate

# 6b) Mapping in minimap2/samtools env
conda activate /mount/britton/Jose/conda-envs/mamba/envs/polypolish_polyca_bwa_bwa_mem2_bowtie2_samtools_pilon_minimap2
if [ ! -s "$OUT/patch_ont.bam.bai" ]; then
  minimap2 -ax map-ont -t "$T" "$FINAL_ASM" "$ONT" \
    | samtools sort -@ "$T" -o "$OUT/patch_ont.bam" -                  
  samtools index "$OUT/patch_ont.bam"
fi
conda deactivate

echo -e "ending:   validate patched  [$(ts)]\n---------------------------"



