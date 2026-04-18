#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --output=/mount/britton/Jose/Jobs/Jason/2025-12-15-In-silico-filtering-last-effort-in-closing-genome-bash-run/workflow/logs/assembly_insilico_filtering_%j.out 
#SBATCH --error=/mount/britton/Jose/Jobs/Jason/2025-12-15-In-silico-filtering-last-effort-in-closing-genome-bash-run/workflow/logs/assembly_insilico_filtering_%j.err  
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=main
#SBATCH --mail-user=jose.cantu@bcm.edu
#SBATCH --mail-type=FAIL,END 

# Load the Anaconda module
module load anaconda3/2023.03-1

# Initialize conda for the shell
eval "$(conda shell.bash hook)" 

# Defining the output directory for easier reference downstream 
OUT=../workflow-output

# Creating output directories for each step besides circulator and bakta otherwise they will not run  
mkdir -p "$OUT"/{reads,flye,map,pilon,quast} 

ts(){ date -Iseconds; }

# 0) Porchop was used here infers ONT adapters directly from reads then trim {May switch out with original Porechop which uses a databases instead of a kmer approach  
conda activate /mount/britton/Jose/conda-envs/porechop_abi 
echo -e "\nstarting: porechop_abi  [$(ts)]\n---------------------------------"
if [ ! -s "$OUT/reads/Pseudo.ont.trim.fastq.gz" ]; then 
	porechop_abi --ab_initio -i ../data/Pseudo_nanopore.fastq.gz -o "$OUT/reads/Pseudo.ont.trim.fastq.gz" -t ${SLURM_CPUS_PER_TASK:-8} 
fi 
conda deactivate
echo -e "ending:   porechop_abi  [$(ts)]\n---------------------------"

# 1) Long-read assembly (Flye)
conda activate /mount/britton/Jose/conda-envs/mamba/envs/autocycler_env
echo -e "\nstarting: flye  [$(ts)]\n---------------------------------"
if [ ! -s "$OUT/flye/assembly.fasta" ]; then 
	flye --nano-raw "$OUT/reads/Pseudo.ont.trim.fastq.gz" --out-dir "$OUT/flye" --genome-size 2.6m --asm-coverage 100 -t ${SLURM_CPUS_PER_TASK:-8} 
fi 
conda deactivate
echo -e "ending:   flye [$(ts)]\n---------------------------------" 

conda activate /mount/britton/Jose/conda-envs/mamba/envs/polypolish_polyca_bwa_bwa_mem2_bowtie2_samtools_pilon_minimap2  
# 2) Map Illumina to the Flye assembly (for polish + coverage filter)
echo -e "\nstarting: bowtie2/samtools  [$(ts)]\n---------------------------------" 
[ -s "$OUT/map/flye.1.bt2" ] || bowtie2-build "$OUT/flye/assembly.fasta" "$OUT/map/flye"
if [ ! -s "$OUT/map/sr.bam.bai" ]; then 
	bowtie2 --very-sensitive -x "$OUT/map/flye" -1 ../data/Pseudo_S670_R1_001.fastq.gz -2 ../data/Pseudo_S670_R2_001.fastq.gz -p ${SLURM_CPUS_PER_TASK:-8} | samtools sort -@ ${SLURM_CPUS_PER_TASK:-8} -o "$OUT/map/sr.bam" -
  	samtools index "$OUT/map/sr.bam"
fi 
echo -e "ending:   bowtie2/samtools  [$(ts)]\n---------------------------"

# 3) Remove low-coverage contigs (≤15× mean short-read depth)
echo -e "\nstarting: coverage filter  [$(ts)]\n---------------------------------"
if [ ! -s "$OUT/pilon/assembly.clean.fasta" ]; then 
	samtools idxstats "$OUT/map/sr.bam" \
 	| awk '$3>0{print $1}' > "$OUT/map/contigs.list"
# compute per-contig coverage and filter
	samtools depth -a "$OUT/map/sr.bam" \
 	| awk '{cov[$1]+=$3; len[$1]++} END{for (c in cov) if (len[c]>0) printf "%s\t%.2f\n", c, cov[c]/len[c]}' \
 	| awk '$2>15{print $1}' > "$OUT/map/keep.contigs" 
	seqtk subseq "$OUT/flye/assembly.fasta" "$OUT/map/keep.contigs" > "$OUT/pilon/assembly.clean.fasta" 
fi 
echo -e "ending:   coverage filter  [$(ts)]\n---------------------------"

# 4) Short-read polish (Pilon; 2–3 rounds if needed)
echo -e "\nstarting: pilon×2  [$(ts)]\n---------------------------------"
export JAVA_TOOL_OPTIONS="-Xmx50g" 
if [ ! -s "$OUT/pilon/polish.done" ]; then 
  for r in 1 2; do
   bowtie2-build "$OUT/pilon/assembly.clean.fasta" "$OUT/pilon/pilon_idx" 
   bowtie2 --very-sensitive -x "$OUT/pilon/pilon_idx" -1 ../data/Pseudo_S670_R1_001.fastq.gz -2 ../data/Pseudo_S670_R2_001.fastq.gz -p ${SLURM_CPUS_PER_TASK:-8} | samtools view -u -F 2308 -s 42.25 | samtools sort -@ ${SLURM_CPUS_PER_TASK:-8} -o "$OUT/pilon/sr.bam" - 
   samtools index "$OUT/pilon/sr.bam" 
   pilon --genome "$OUT/pilon/assembly.clean.fasta" --frags "$OUT/pilon/sr.bam" --threads ${SLURM_CPUS_PER_TASK:-8} --outdir "$OUT/pilon" --output "pilon${r}"
   test -s "$OUT/pilon/pilon${r}.fasta" || { echo "FATAL pilon round ${r} failed"; exit 2; }
   mv "$OUT/pilon/pilon${r}.fasta" "$OUT/pilon/assembly.clean.fasta" 
  done
  : > "$OUT/pilon/polish.done" # mark success
fi 
conda deactivate 
echo -e "ending:   pilon×2  [$(ts)]\n---------------------------"

echo -e "\nstarting: circlator  [$(ts)]\n---------------------------------"
conda activate /mount/britton/Jose/conda-envs/mamba/envs/circlator_env 
# 5) Circularization attempt
if [ ! -s "$OUT/circulator/06.fixstart.fasta" ]; then 
	circlator all --threads ${SLURM_CPUS_PER_TASK:-8} "$OUT/pilon/assembly.clean.fasta" "$OUT/reads/Pseudo.ont.trim.fastq.gz" "$OUT/circulator" 
fi 
# If circularized, use Circlator’s output; otherwise keep assembly.clean.fasta
conda deactivate 
echo -e "ending:   circlator  [$(ts)]\n---------------------------"

# 6) QC and annotation
conda activate /mount/britton/Jose/conda-envs/quast
FINAL_ASM="$OUT/pilon/assembly.clean.fasta"; [ -s "$OUT/circulator/06.fixstart.fasta" ] && FINAL_ASM="$OUT/circulator/06.fixstart.fasta"
echo -e "\nstarting: quast  [$(ts)]\n---------------------------------"
[ -s "$OUT/quast/report.txt" ] || quast -o "$OUT/quast" -t ${SLURM_CPUS_PER_TASK:-8} "$FINAL_ASM"
conda deactivate 
echo -e "ending:   quast  [$(ts)]\n---------------------------"


conda activate /mount/britton/Jose/conda-envs/bakta_env
echo -e "\nstarting: bakta  [$(ts)]\n---------------------------------"
export BAKTA_DB_PATH="/mount/britton/leap/dbs/bakta/db" 
[ -s "$OUT/bakta/clean_assembly.gff3" ] || bakta --db "$BAKTA_DB_PATH" --prefix clean_assembly --threads ${SLURM_CPUS_PER_TASK:-8} --output "$OUT/bakta" "$FINAL_ASM"
conda deactivate 
echo -e "ending:   bakta  [$(ts)]\n---------------------------"
