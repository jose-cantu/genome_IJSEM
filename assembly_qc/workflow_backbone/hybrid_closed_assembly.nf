#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

  main:

  // Input CSV: sample,r1,r2,ont
  def csvPath = params.index?.trim()
  def csvFile = file(csvPath)
  if( !csvFile.exists() )
      error "CSV not found on this node: ${csvFile}"

  // Create a single-value channel carrying the CSV file, then split rows
  reads = Channel.value(csvFile)
      .splitCsv(header: true)
      // drop bad/incomplete lines early
      .filter { row -> row.sample && row.r1 && row.r2 && row.ont }
      .map { row ->
          tuple(
              row.sample as String,
              file(row.r1),
              file(row.r2),
              file(row.ont)
          )
      }

// Debug dump (non-destructive)
reads.view { "READS: $it" }
 
  // QC and Trimming 
  trimmed  = FASTP(reads)
  ont_trim = PORECHOP_ABI(trimmed)
  ont_qc   = NANOPLOT(ont_trim)
  
  // Autocycler assembly -> rotation 
  asm_lr   = AUTOCYCLER(ont_qc)
  rot      = DNAAPLER_ROTATE(asm_lr)

  // Medaka A/B test to see if Medaka improves overall quality of final genome assembly or not (recent papers looked at down vote others mixed signals) 
  // only execute if params.medaka_probe is true 
  def yak_ab = nextflow.Channel.empty()  
  if (params.medaka_probe) { 
     // Apply Medaka to the rotated Autocycler assembly 
     medaka_polish = MEDAKA_PROBE(rot)
     // compare with/without Medka using Yak 
     yak_ab = YAK_AB_TEST(
         medaka_polish.map { sample, r1t, r2t, ont_t, fastp_html, fastp_json, nanoplot_dir, asm_no_medaka, asm_medaka ->
            tuple(sample, r1t, r2t, asm_no_medaka, asm_medaka)
         }
     )
  }

  // Coverage on primary alignments; produce coverage.tsv and mean_depth.txt 
  cov      = COV_PRIMARY(rot)

  // Parse mincov.txt -> numeric for depth aware triage; keep tuple shape + mincov at end and ignore coveerage.tsv and mean_deapth.txt 
  cov_num = cov.map { sample, r1t, r2t, ont_t, fastp_html, fastp_json, nanoplot_dir, asm_fa, asm_gfa, coverage_tsv, mean_depth_txt, mincov_txt ->
    def v = (mincov_txt.text?.trim() ?: '0').toDouble()
    tuple(sample, r1t, r2t, ont_t, fastp_html, fastp_json, nanoplot_dir, asm_fa, asm_gfa, v)
  }

  // Depth-aware branching (add a catch-all)
  br = cov_num.branch { sample, r1t, r2t, ont_t, fastp_html, fastp_json, nanoplot_dir, asm_fa, asm_gfa, mincov ->
    low  : (mincov < 5)
    mid  : (mincov >= 5 && mincov < 25)
    high : (mincov >= 25)
    other: true    // safe bucket per docs
  }
  
  // create one empty channel as a fallback
  def emptyChan = nextflow.Channel.empty()
  
  // assign the branch channels, falling back to emptyChan
  def lowChan  = br.low  ?: emptyChan
  def midChan  = br.mid  ?: emptyChan
  def highChan = br.high ?: emptyChan

  // call polishing only when there is data; otherwise return the empty channel
  def pol_low  = POLYPOLISH_CAREFUL_ONLY(lowChan) 
  def pol_mid  = POLY_CAREFUL_PLUS_PYPOLCA(midChan) 
  def pol_high = POLY_DEFAULT_PLUS_PYPOLCA(highChan) 

  final_asm = pol_low.mix(pol_mid).mix(pol_high)
  
  // QC and annotation 
  quast  = QUAST(final_asm)
  checkm1 = CHECKM1(quast)
  busco  = BUSCO(checkm1)
  yak    = YAK(busco)
  bakta  = BAKTA(yak)

  // Extract 16S rRNA sequences for Sanger confirmation (Barrnap)
  barr = BARRNAP_16S(bakta)

 
  // GTDB-Tk classification and phylogenomics (de novo MSA and IQ-TREE) 
  gtdb   = GTDBTK(bakta)
  gtdb_tree = GTDB_CLASSIFY_TREE(gtdb)

  // Build refs.txt from GTDB neighbors
  // Build the FastANI reference panel from the GTDB classification summary is what I'm doing.... 
  panel = BUILD_FASTANI_PANEL(
    gtdb_tree.map { sample, asm_fa, asm_gfa, gtdbtk_dir,
                      quast_dir, checkm1_dir, busco_dir, yak_qv,
                      bakta_dir, class_tree, class_summary ->
    tuple(sample, class_summary)
      }
   )

  // Prepare (sample, asm_final_fa) pairs from the barr channel
  asm_for_fastani = barr.map { sample, asm_fa, asm_gfa, quast_dir,
                                 checkm1_dir, busco_dir, yak_qv, bakta_dir, f16S ->
      tuple(sample, asm_fa)
    }
   
  fastani_input = asm_for_fastani.join(panel).map { s,a,refs -> tuple(s, a, refs) }

  // enforce that FastANI must have references
  if (params.fastani_required) {
      fastani_input = fastani_input.filter { s,a,refs -> refs.exists() && refs.size() > 0 }
      if (!fastani_input) error "FastANI required but no reference genomes were obtained."
    }

  // Run FastANI
  fastani = FASTANI_BATCH(fastani_input)

  // Prepare dDDH bundle (assemblies + annotation + FastANI)
  tygs_input = barr.join(fastani).map { sample, asm_fa, asm_gfa, quast_dir,
                                      checkm1_dir, busco_dir, yak_qv,
                                      bakta_dir, f16S, fastani_tsv ->
  tuple(sample, asm_fa, bakta_dir, quast_dir, checkm1_dir,
        busco_dir, yak_qv, f16S, fastani_tsv)
  }
  tygs = TYGS_EXPORT(tygs_input)
  
  // Jasons tree workflow gtdb -> gtdbtree -> itol 
  denovo_msa = GTDB_DENOVO_ALIGN(gtdb)
  // publication denovo grade tree gtdb -> gtdb denovo align -> iqtree2
  denovo_brief = denovo_msa.map { sample, asm_fa, asm_gfa, gtdbtk_dir, quast_dir,
                                  checkm1_dir, busco_dir, yak_qv, bakta_dir, denovo_dir, msa_gz, ftree ->
	tuple(sample, msa_gz, ftree) 
	} 
  iqtree = IQTREE2(denovo_brief)
  itol = ITOL(gtdb_tree) 
  band   = BANDAGE_IMG(gtdb)
  MULTIQC(band)
}

/* ------------------- FASTP Process ---------------------- */

// Illumina QC
process FASTP {
  conda '/mount/britton/Jose/conda-envs/fastp'
  tag { sample }
  publishDir { "${params.outdir}/fastp/${sample}/" }, mode: 'copy'

  input:
  tuple val(sample), path(r1), path(r2), path(ont)

  output:
  tuple val(sample),
        path("${sample}.R1.fastp.fq.gz"),
        path("${sample}.R2.fastp.fq.gz"),
        path(ont),
        path("${sample}.fastp.html"),
        path("${sample}.fastp.json")

  script:
  """
  # fastp: trim adapters/low-quality bases; emit HTML+JSON QC reports
  fastp \
    -i "${r1}" -I "${r2}" \
    -o "${sample}.R1.fastp.fq.gz" \
    -O "${sample}.R2.fastp.fq.gz" \
    -h "${sample}.fastp.html" \
    -j "${sample}.fastp.json" \
    -w ${task.cpus}
  """
}

/* --------------- PORECHOP_ABI Process ---------------- */

// ONT adapter trim (ab initio)
process PORECHOP_ABI {
  conda '/mount/britton/Jose/conda-envs/porechop_abi'
  tag { sample }
  publishDir { "${params.outdir}/porechop_abi/${sample}" }, mode: 'copy'

  input:
  tuple val(sample), path(r1t), path(r2t), path(ont),
        path(fastp_html), path(fastp_json)

  output:
  tuple val(sample),
        path(r1t), path(r2t),
        path("${sample}.ont.trim.fastq.gz"),
        path(fastp_html), path(fastp_json)

  script:
  """
  # Porechop_ABI: infer ONT adapters directly from reads, then trim
  porechop_abi -abi \
    -i "${ont}" \
    -o "${sample}.ont.trim.fastq.gz" \
    -t ${task.cpus}
  """
}

// 3) NanoPlot
process NANOPLOT {
  conda '/mount/britton/Jose/conda-envs/NanoPlot'
  tag { sample }
  publishDir { "${params.outdir}/NanoPlot/${sample}/" }, mode: 'copy'

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json),
        path("nanoplot")

  script:
  """
  # NanoPlot: QC plots for ONT read length/quality
  mkdir -p nanoplot
  NanoPlot \
    --fastq "${ont_trim}" \
    --threads ${task.cpus} \
    -o nanoplot
  """
}

// 4) Autocycler (long read first assembly approach; Using 9 different assemblers on four subsets, compressing, clustering, trimming, and resolving and combining contigs and connecting the best contigs from each to create a unified circucle final assembly and then dnaapler rotates circular replicons to dnaA/repA at next step)
process AUTOCYCLER {
  tag { sample }
  publishDir { "${params.outdir}/Autocycler/${sample}/" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/mamba/envs/autocycler_env"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path("assembly.fasta"), path("assembly.gfa")

  shell:
  '''
  autocycler --version || true 2>&1 | tee -a versions.log
  parallel    --version || true 2>&1 | tee -a versions.log

  READS="!{ont_trim}"
  THREADS=!{params.autocycler_threads}
  JOBS=!{params.autocycler_jobs}
  READ_TYPE="!{params.autocycler_readtype}"
  MAX_TIME="!{params.autocycler_timeout}"
  MINDEP="!{params.autocycler_mindepth}"



  # 1) Estimate genome size
  genome_size=$(autocycler helper genome_size --reads "$READS" --threads "$THREADS")

  # 2) Subsample long reads into four non-overlapping subsets
  autocycler subsample --reads "$READS" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

  # 3) Build job list for nine assemblers across 4 subsets
  mkdir -p assemblies
  : > assemblies/jobs.txt
  for assembler in raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu; do
    for i in 01 02 03 04; do
      echo "autocycler helper $assembler --reads subsampled_reads/sample_${i}.fastq --out_prefix assemblies/${assembler}_${i} --threads $THREADS --genome_size $genome_size --read_type $READ_TYPE --min_depth_rel $MINDEP" >> assemblies/jobs.txt
    done
  done

  # 4) Run helper jobs in parallel with timeout
  set +e
  nice -n 19 parallel --jobs "$JOBS" --joblog assemblies/joblog.tsv --results assemblies/logs --timeout "$MAX_TIME" < assemblies/jobs.txt
  set -e

  # Optional: free space
  rm -f subsampled_reads/*.fastq || true

  # 5–7) Compress, cluster, trim & resolve
  autocycler compress -i assemblies -a autocycler_out
  autocycler cluster  -a autocycler_out
  for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim    -c "$c"
    autocycler resolve -c "$c"
  done

  # 8) Combine into final assembly
  autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

  # Normalize outputs to declared names
  cp -f autocycler_out/consensus_assembly.fasta assembly.fasta
  if [[ -f autocycler_out/consensus_assembly.gfa ]]; then
    cp -f autocycler_out/consensus_assembly.gfa assembly.gfa
  else
    : > assembly.gfa
  fi
  '''
  }

// 5) Rotation of Assembled Genome using DNAAPLER to proper orientation 
process DNAAPLER_ROTATE { 
  tag { sample }
  publishDir { "${params.outdir}/Dnaapler/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/mamba/envs/dnaapler_env"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path("assembly.fasta"), path(asm_gfa)

  script:
  """
  dnaapler --version || true 2>&1 | tee -a versions.log

  # Rotate to dnaA/repA if found; otherwise pass through unchanged

  cp -f "${asm_fa}" in.fasta
  dnaapler rotate \
    --fasta in.fasta \
    --out rotated.fasta \
    --genes "${params.dnaapler_genes ?: 'dnaA,repA' }" \
    --threads ${task.cpus} || cp -f in.fasta rotated.fasta

  cp -f rotated.fasta assembly.fasta
  """
}

// Testing to see if MEDAKA Improves overall quality of GENOME ASSEMBLY mixed signals from latest papers testing 
process MEDAKA_PROBE {
  tag { sample }
  publishDir { "${params.outdir}/medaka/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/medaka"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_no_medaka), path(asm_gfa)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_no_medaka),  // for A/B reference
        path("assembly.medaka.fasta")

  script:
  """
  medaka_consensus --version || true 2>&1 | tee -a versions.log

  model="${ params.medaka_model ?: 'r104_e81_sup_g610' }"  # set to your R10.4 sup model
  # medaka_consensus can take raw reads + draft
  medaka_consensus \
    -i "${ont_trim}" \
    -d "${asm_no_medaka}" \
    -o medaka_out \
    -t ${task.cpus} \
    -m "${model}"

  # Produce FASTA
  if [[ -f medaka_out/consensus.fasta ]]; then
    cp -f medaka_out/consensus.fasta "assembly.medaka.fasta"
  elif [[ -f medaka_out/consensus.fna ]]; then
    cp -f medaka_out/consensus.fna "assembly.medaka.fasta"
  else
    # fall back: pass original assembly if medaka did not produce consensus
    cp -f "${asm_no_medaka}" "assembly.medaka.fasta"
  fi
  """
}

// Yak Testing seeing if MEDAKA Iproves or not comparing with vs without 
process YAK_AB_TEST {
  tag { sample }
  publishDir { "${params.outdir}/yak_ab/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/yak"

  input:
  tuple val(sample), path(r1t), path(r2t), path(asm_no_medaka), path(asm_medaka)

  output:
  path("yak_ab.tsv"), path("yak_no_medaka.txt"), path("yak_medaka.txt")

  script:
  """
  yak --version || true 2>&1 | tee -a versions.log

  yak count -o ${sample}.yak "${r1t}" "${r2t}"

  yak qv -t ${task.cpus} ${sample}.yak "${asm_no_medaka}" > yak_no_medaka.txt
  ADJ_A=\$(awk '\$1=="QV"{print \$3}' yak_no_medaka.txt)

  if [[ -s "${asm_medaka}" ]]; then
    yak qv -t ${task.cpus} ${sample}.yak "${asm_medaka}" > yak_medaka.txt
    ADJ_B=\$(awk '\$1=="QV"{print \$3}' yak_medaka.txt)
  else
    : > yak_medaka.txt
    ADJ_B="NA"
  fi

  DELTA=\$(awk -v a="${ADJ_A}" -v b="${ADJ_B}" 'BEGIN{ if(b=="NA"||a==""){print "NA"; exit} printf("%.2f", b-a) }')
  REC=\$(awk -v d="${DELTA}" 'BEGIN{ if(d=="NA") print "no_medaka_only"; else if(d+0>=0.5) print "use_medaka"; else print "skip_medaka" }')

  printf "sample\tadj_qv_no_medaka\tadj_qv_medaka\tdelta\trecommendation\n" > yak_ab.tsv
  printf "%s\t%s\t%s\t%s\t%s\n" "${sample}" "${ADJ_A}" "${ADJ_B}" "${DELTA}" "${REC}" >> yak_ab.tsv
  """
}

// 5) Coverage (primary alignments) -> mincov.txt
process COV_PRIMARY {
  tag { sample }
  publishDir { "${params.outdir}/Cov_Primary/${sample}/" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/polypolish_polyca_bwa_bwa_mem2_samtools"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa),
        path("coverage.tsv"), path("mean_depth.txt"), path("mincov.txt") 

  script:
  """
  bwa-mem2 index "${asm_fa}"
  bwa-mem2 mem -t "${task.cpus}" "${asm_fa}" "${r1t}" "${r2t}" \
    | samtools view -F 2304 -b - \
    | samtools sort -o sr.primary.bam

  samtools index sr.primary.bam
  samtools coverage -o coverage.tsv sr.primary.bam

  awk 'NR>1{len=\$3-\$2+1; sum+=len*\$7; L+=len} END{if(L>0) printf("%.2f\\n", sum/L); else print "0.00"}' coverage.tsv > mean_depth.txt
  awk 'NR>1{if(min==""||\$7<min)min=\$7} END{printf("%.2f\\n", (min=="")?0:min)}' coverage.tsv > mincov.txt
  """

}

// 6a) <5× : Polypolish --careful only
process POLYPOLISH_CAREFUL_ONLY {
  tag { sample }
  publishDir { "${params.outdir}/polishing_low_depth/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/polypolish_polyca_bwa_bwa_mem2_samtools"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa), val(mincov)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path("assembly.polished.fasta"), path(asm_gfa)

  script:
  """
  # Build index for all-alignments mapping required by Polypolish
  bwa-mem2 index "${asm_fa}"

  # Map each short-read end with -a (emit multi-mapping) to produce all-alignments SAMs
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r1t}" > R1.all.sam
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r2t}" > R2.all.sam

  # Polypolish: short-read consensus polishing (conservative mode)
  polypolish polish --careful "${asm_fa}" R1.all.sam R2.all.sam > "assembly.polished.fasta"
  """
}

// 6b) 5–25× : Polypolish --careful -> Pypolca --careful
process POLY_CAREFUL_PLUS_PYPOLCA {
  tag { sample }
  publishDir { "${params.outdir}/polishing_mid_depth/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/polypolish_pypolca_bwa-mem2_samtools"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa), val(mincov)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path("assembly.polished.fasta"), path(asm_gfa)

  script:
  """
  # Index assembly and create all-alignments SAMs for Polypolish
  bwa-mem2 index "${asm_fa}"
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r1t}" > R1.all.sam
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r2t}" > R2.all.sam

  # First pass: Polypolish in conservative mode
  polypolish polish --careful "${asm_fa}" R1.all.sam R2.all.sam > pp1.fasta

  # Second pass: pypolca with --careful thresholds to suppress low-depth false positives
  pypolca run \
    -a pp1.fasta \
    -1 "${r1t}" -2 "${r2t}" \
    -t ${task.cpus} \
    --careful \
    -o polca_out
  cp polca_out/pypolca_corrected.fasta "assembly.polished.fasta"
  """
}

// 6c) >25× : Polypolish (default) -> Pypolca --careful
process POLY_DEFAULT_PLUS_PYPOLCA {
  tag { sample }
  publishDir { "${params.outdir}/polishing_high_depth/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/polypolish_polyca_bwa_bwa_mem2_samtools"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_fa), path(asm_gfa), val(mincov)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path("assembly.polished.fasta"), path(asm_gfa)

  script:
  """
  # Index assembly and generate Polypolish all-alignments SAMs
  bwa-mem2 index "${asm_fa}"
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r1t}" > R1.all.sam
  bwa-mem2 mem -a -t ${task.cpus} "${asm_fa}" "${r2t}" > R2.all.sam

  # First pass: Polypolish in default mode (suitable at ≥25× depth)
  polypolish polish "${asm_fa}" R1.all.sam R2.all.sam > pp1.fasta

  # Second pass: pypolca --careful (robust variant for residual SNVs/indels)
  pypolca run \
    -a pp1.fasta \
    -1 "${r1t}" -2 "${r2t}" \
    -t ${task.cpus} \
    --careful \
    -o polca_out
  cp polca_out/pypolca_corrected.fasta "assembly.polished.fasta"
  """
}

// 7) QUAST
process QUAST {
  tag { sample }
  publishDir { "${params.outdir}/QUAST/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/quast"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa),
        path("quast")

  script:
  """
  # QUAST: assembly structural metrics
  quast \
    -o quast \
    -t ${task.cpus} \
    "${asm_final_fa}"
  """
}

// 8) CheckM v1.2.2
process CHECKM1 {
  tag { sample }
  publishDir { "${params.outdir}/Checkm1/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/checkm"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa), path(quast_dir)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path("checkm1")

  script:
  """
  # Prepare genome set for CheckM (expects a directory of FASTA files)
  mkdir -p genomes
  cp "${asm_final_fa}" "genomes/${sample}.fasta"

  # CheckM v1 lineage workflow: completeness/contamination
  checkm lineage_wf \
    -x fasta \
    -t ${task.cpus} \
    genomes checkm1
  """
}

// 9) BUSCO (odb10)
process BUSCO {
  tag { sample }
  publishDir { "${params.outdir}/busco/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/busco"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir)

  output:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path("busco")

  script:
  """
  # BUSCO: single-copy ortholog completeness (offline)
  busco \
    -i "${asm_final_fa}" \
    -m genome \
    -l "${params.busco_lineage}" \
    -o busco \
    -c ${task.cpus} \
    --offline
  """
}

// 10) Yak QV
process YAK {
  tag { sample }
  publishDir { "${params.outdir}/yak/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/yak"

  input:
  tuple val(sample),
        path(r1t), path(r2t), path(ont_trim),
        path(fastp_html), path(fastp_json), path(nanoplot_dir),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir)

  output:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir),
        path("yak.qv.txt"),
        path(fastp_html), path(fastp_json), path(nanoplot_dir)

  script:
  """
  # Build k-mer DB from short reads (yak count)
  yak count \
    -o "${sample}.yak" \
    "${r1t}" "${r2t}"

  # Estimate assembly QV by k-mer concordance (yak qv)
  yak qv \
    -t ${task.cpus} \
    "${sample}.yak" \
    "${asm_final_fa}" \
    > "yak.qv.txt"
  """
}

// 11) Bakta
process BAKTA {
  tag { sample }
  publishDir { "${params.outdir}/bakta/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/bakta_env"

  input:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(fastp_html), path(fastp_json), path(nanoplot_dir)

  output:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path("bakta")

  script:
  """
  # Bakta: standardized prokaryotic genome annotation
  bakta \
    --db "${params.bakta_db}" \
    --threads ${task.cpus} \
    --output bakta \
    --prefix "${sample}" \
    "${asm_final_fa}"
  """
}


process BARRNAP_16S {
  conda '/mount/britton/Jose/conda-envs/mamba/envs/barrnap_env'
  tag { sample }
  publishDir { "${params.outdir}/16S/${sample}" }, mode: 'copy'

  input:
  tuple val(sample), path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir)

  output:
  tuple val(sample), path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir), path('16S.fna')

  script:
  """
  barrnap --version || true 2>&1 | tee -a versions.log
  samtools --version || true 2>&1 | tee -a versions.log

  barrnap --kingdom bac --threads ${task.cpus} "${asm_final_fa}" > rRNA.gff
  samtools faidx "${asm_final_fa}"
  awk '\$3=="rRNA" && \$9~"16S"{print \$1":"\$4"-"\$5}' rRNA.gff > coords.txt
  if [[ -s coords.txt ]]; then
    xargs -a coords.txt -I{} samtools faidx "${asm_final_fa}" {} > 16S.fna
  else
    : > 16S.fna
  fi
  """ 
  }


// 12a) Build FastANI reference panel from GTDB summary
process BUILD_FASTANI_PANEL {
  tag { sample }
  publishDir { "${params.outdir}/fastani_panel/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/ncbi_datasets"

  input:
    tuple val(sample), path(class_summary)

  output:
    tuple val(sample), path("refs.txt")

  // NOTE: use triple double quotes + Groovy defaults inside !{ }
  shell:
  '''
  set -euo pipefail
  
  # 1) Extract candidate accessions by header name (no fixed column numbers)
  awk -F '\t' '
    function emit_token(tok,    p,acc) {
      p = index(tok,"GCF_"); if (!p) p = index(tok,"GCA_");
      if (p) {
        acc = substr(tok,p)
              # trim at first space/comma/semicolon (no regex)
        sp = index(acc," "); cm = index(acc,","); sc = index(acc,";"); cut = 0
        if (sp && (!cut || sp<cut)) cut = sp
        if (cm && (!cut || cm<cut)) cut = cm
        if (sc && (!cut || sc<cut)) cut = sc
        if (cut) acc = substr(acc,1,cut-1)
        if (index(acc,".")>0) print acc
      }
    }
    NR==1{
      for (i=1;i<=NF;i++) h[\$i]=i

      pr = 0
      if (h["fastani_reference"])                 pr = h["fastani_reference"]
      else if (h["closest_placement_reference"])  pr = h["closest_placement_reference"]
      else if (h["closest_reference"])            pr = h["closest_reference"]
      else if (h["closest_genome"])               pr = h["closest_genome"]

      orr = 0
      if (h["other_related_references(genome_id,species_name,radius,ANI,AF)"])
           orr = h["other_related_references(genome_id,species_name,radius,ANI,AF)"]
      else if (h["other_related_references"])
           orr = h["other_related_references"]

      next
    }
    {
      if (pr  && \$pr!=""  && \$pr!="NA")  emit_token(\$pr)
      if (orr && \$orr!="" && \$orr!="NA") {
          s = \$orr; gsub(",", " ", s); gsub(";", " ", s)
          n = split(s,a," ")
          for (i=1;i<=n;i++) if (a[i]!="") emit_token(a[i])
      }
    }
  ' "!{class_summary}" | sort -u | head -n !{ params.fastani_max_refs ?: 20 } > accessions.txt
  mkdir -p panel meta

  # 2) Download each accession and extract a FASTA
  while read -r acc; do
    datasets summary genome accession "\$acc" --as-json-lines > "meta/\${acc}.jsonl" || true
    datasets download genome accession "\$acc" --include genome --filename "\${acc}.zip" || continue
    unzip -oq "\${acc}.zip" -d "\${acc}" || true
    f=\$(find "\${acc}/ncbi_dataset/data" -type f \\( -name "*_genomic.fna" -o -name "*.fna" \\) | head -n1 || true)
    if [[ -n "\$f" ]]; then
      cp -f "\$f" "panel/\${acc}.fna"
    fi
    rm -rf "\${acc}.zip" "\${acc}"
  done < accessions.txt

  # 3) Prefer type-material assemblies if any are available
  if [[ "!{ params.fastani_prefer_type ?: true }" == "true" ]]; then
    : > refs.type.txt
    while read -r acc; do
      if jq -e -r 'select(.assembly_info.relation_to_type_material // "" | test("type"; "i"))' "meta/\${acc}.jsonl" >/dev/null 2>&1; then
        [[ -s "panel/\${acc}.fna" ]] && echo "panel/\${acc}.fna" >> refs.type.txt
      fi
    done < accessions.txt
    if [[ -s refs.type.txt ]]; then
      sort -u refs.type.txt > refs.txt
      exit 0
    fi
  fi

  # 4) Otherwise, use all downloaded neighbours
  if compgen -G "panel/*.fna" > /dev/null; then
    realpath panel/*.fna > refs.txt
  else
    : > refs.txt
  fi
  '''
}


process FASTANI_BATCH {
  conda '/mount/britton/Jose/conda-envs/fastani'
  tag { sample }
  publishDir { "${params.outdir}/fastani/${sample}" }, mode: 'copy'

  input:
  tuple val(sample), path(asm_final_fa), path(ref_list)

  output:
  tuple val(sample), path('fastani.tsv')

  script:
  """
  { fastANI --version || true; } 2>&1 | tee -a versions.log

  if [[ ! -s "${ref_list}" ]]; then
    echo "WARN: empty ref list; writing minimal fastani.tsv" | tee -a versions.log
    printf "query\\treference\\tani\\tfract\\tqlen\\trlen\\taln_frags\\n" > fastani.tsv
    exit 0
  fi


  fastANI --threads ${task.cpus} \
          --query "${asm_final_fa}" \
          --refList "${ref_list}" \
          --output fastani.tsv
  """
}



// Export TYGS needs ..... 
process TYGS_EXPORT {
  tag { sample }
  publishDir { "${params.outdir}/tygs/${sample}" }, mode: 'copy'

  input:
  tuple val(sample), path(asm_final_fa), path(bakta_dir),
        path(quast_dir), path(checkm1_dir), path(busco_dir),
        path(yak_qv), path(f16S), path(fastani_tsv)

  output:
  path("${sample}.tygs_bundle.tgz")

  script:
  """
  set -euo pipefail

  tar -czf ${sample}.tygs_bundle.tgz \
    -C "\$(dirname \"${asm_final_fa}\")" "\$(basename \"${asm_final_fa}\")" \
    -C "${bakta_dir}" . \
    -C "${quast_dir}" . \
    -C "${checkm1_dir}" . \
    -C "${busco_dir}" . \
    -C "\$(dirname \"${yak_qv}\")" "\$(basename \"${yak_qv}\")" \
    -C "\$(dirname \"${f16S}\")" "\$(basename \"${f16S}\")" \
    -C "\$(dirname \"${fastani_tsv}\")" "\$(basename \"${fastani_tsv}\")" || true
  """
}

// 12) GTDB-Tk 2.1.1 (R207_v2)
process GTDBTK {
  tag { sample }
  publishDir { "${params.outdir}/gtdbtk/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/gtdbtk"

  input:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir)

  output:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa), path("gtdbtk"),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir)

  script:
  """
  # Prepare genome list for GTDB-Tk classify_wf
  mkdir -p genomes
  cp "${asm_final_fa}" "genomes/${sample}.fasta"

  # Point GTDB-Tk to the desired reference release
  export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

  # GTDB-Tk: assign taxonomy using GTDB reference
  gtdbtk classify_wf \
    --genome_dir genomes \
    --out_dir gtdbtk \
    -x fasta \
    --cpus ${task.cpus}
  """
}

// 13a) classifying tree in gtdb file harvesting step to extract gtdb classification tree and summary for downstream process 
process GTDB_CLASSIFY_TREE {
  tag { sample } 
  publishDir { "${params.outdir}/trees/gtdb_classify/${sample}" }, mode: 'copy' 
  conda "/mount/britton/Jose/conda-envs/gtdbtk"

  input:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir), 
	path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	path(bakta_dir)

  output:
  tuple val(sample), 
  	path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	path(bakta_dir), path("gtdb_classify.tree"), path("gtdb_classify.summary.tsv")

  script:
  """
  echo "== Listing gtdbtk_dir ==" >&2
  find "${gtdbtk_dir}" -maxdepth 5 -type f -printf '%P\\n' >&2 || true

  # Prefer typical GTDB-Tk summary names, fall back to any *summary.tsv
  S=\$(find "${gtdbtk_dir}" -maxdepth 5 -type f -name 'gtdbtk.*summary.tsv' -print -quit 2>/dev/null || true)
  if [[ -z "\$S" ]]; then
    S=\$(find "${gtdbtk_dir}" -maxdepth 5 -type f -name '*summary.tsv' -print -quit 2>/dev/null || true)
  fi

  # Tree file
  T=\$(find "${gtdbtk_dir}" -maxdepth 5 -type f -name '*.tree' -print -quit 2>/dev/null || true)

  # Copy or create placeholders
  if [[ -n "\$T" ]]; then
    cp -f "\$T" gtdb_classify.tree
  else
    printf "(%s);\n" "${sample}" > gtdb_classify.tree
  fi

  if [[ -n "\$S" ]]; then
    cp -f "\$S" gtdb_classify.summary.tsv
  else
    echo "WARNING: No GTDB-Tk summary found under ${gtdbtk_dir} — writing empty file" >&2
    : > gtdb_classify.summary.tsv
  fi
  """
}

// 13b) gtdb devo alignment tree build de-novo marker MSA (including GTDB refs) 
process GTDB_DENOVO_ALIGN { 
    tag { sample } 
    publishDir { "${params.outdir}/trees/gtdb_denovo/${sample}" }, mode: 'copy' 
    conda "/mount/britton/Jose/conda-envs/gtdbtk" 

    input: 
    tuple val(sample),
    	  path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
	  path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	  path(bakta_dir)

    output: 
    tuple val(sample),
    	  path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
	  path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	  path(bakta_dir), 
	  path("denovo"), // full gtdbtk output directory 
	  path("user_msa.fasta.gz"), // masked concatenated MSA 
	  path("fasttree.tree") // quick tree reference 
    
    shell:
    """ 
    export GTDBTK_DATA_PATH="!{params.gtdbtk_db}" 
    mkdir -p genomes denovo 
    cp -f "!{asm_final_fa}" "genomes/!{sample}.fna"
    abs=\$(realpath "genomes/!{sample}.fna")
    printf "%s\t!{sample}\n" "\$abs" > batchfile.tsv
    gtdbtk de_novo_wf --genome_dir genomes --bacteria --outgroup_taxon p__Chloroflexota --out_dir denovo --cpus !{task.cpus} -x fna 

    M=\$(find denovo -type f -name "*msa*.fasta.gz" -print -quit || true) 
    T=\$(find denovo -type f -name "*.tree" -print -quit || true) 
    [[ -n "\$M" ]] && cp -f "\$M" user_msa.fasta.gz || : > user_msa.fasta.gz 
    [[ -n "\$T" ]] && cp -f "\$T" fasttree.tree || printf "(%s);\n" "!{sample}" > fasttree.tree 
    
    """ 
    
} 

// 13c) reinfer tree w IQTREE2 (ModelFinder + UFBoost +SH-aLRT) 
process IQTREE2 { 
   tag { sample } 
   publishDir { "${params.outdir}/trees/iqtree/${sample}" }, mode: 'copy' 
   conda "/mount/britton/Jose/conda-envs/iqtree" 
   
   input:
   tuple val(sample), 
   	 path(msa_gz), path(ftree) 

   output:
   tuple val(sample),
         path("iqtree/${sample}.treefile"),
	 path("iqtree/${sample}.iqtree.log")

   script:
   """
   mkdir -p iqtree 
   gunzip -c "${msa_gz}" > msa.fasta 
   iqtree2 -s msa.fasta -m MFP -B 1000 --alrt 1000 -T AUTO -pre "iqtree/${sample}" 
   mv "iqtree/${sample}.log" "iqtree/${sample}.iqtree.log" || true 

   """ 
}


// 13d) Use ITOL on Jason's tree to create netwick tree to upload 
process ITOL {
   tag { sample } 
   publishDir { "${params.outdir}/trees/itol/${sample}" }, mode: 'copy'
   conda "/mount/britton/Jose/conda-envs/gtdbtk" 

   input: 
   tuple val(sample),
   	 path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
	 path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	 path(bakta_dir),
	 path(class_tree), path(class_summary)

   output: 
   tuple val(sample),
   	 path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
	 path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
	 path(bakta_dir),
	 path("itol") 

   script:
   """
   mkdir -p itol 
   gtdbtk convert_to_itol \
   	--input_tree "${class_tree}" \
	--output_tree "itol/${sample}.itol.tree" 
   # copy over the summary file into itol direcotyr for reference 
   cp "${class_summary}" "itol/${sample}.summary.tsv" || true
   """ 	
}

// 14) Bandage graph
process BANDAGE_IMG {
  tag { sample }
  publishDir { "${params.outdir}/bandage/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/bandage"

  input:
  tuple val(sample),
        path(asm_final_fa), path(asm_gfa), path(gtdbtk_dir),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir)

  output:
  tuple val(sample),
        path(asm_final_fa), path(gtdbtk_dir),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir),
        path("${sample}.assembly.gfa.png")

  script:
  """
  # Bandage (headless): render assembly graph to PNG
  Bandage image \
    "${asm_gfa}" \
    "${sample}.assembly.gfa.png" \
    --height 1200 --width 1600
  """
}

// 15) MultiQC
process MULTIQC {
  tag { sample }
  publishDir { "${params.outdir}/multiqc/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/multiqc"

  input:
  tuple val(sample),
        path(asm_final_fa), path(gtdbtk_dir),
        path(quast_dir), path(checkm1_dir), path(busco_dir), path(yak_qv),
        path(bakta_dir), path(bandage_png)

  output:
  path("multiqc_report.html")

  script:
  """
  # Stage key tool outputs for aggregation
  mkdir -p in
  mkdir -p in/yak
  mkdir -p in/bandage
  cp -r "${quast_dir}"    in/quast
  cp -r "${busco_dir}"    in/busco
  cp -r "${checkm1_dir}"  in/checkm1 || true
  cp -r "${gtdbtk_dir}"   in/gtdbtk  || true
  cp -f "${asm_final_fa}" "in/${sample}.assembly.fasta" || true

  # Yak QV -> simple TSV
  if [[ -f "${yak_qv}" ]]; then
    cp -f "${yak_qv}" "in/yak/${sample}.yak.qv.txt"
    QV=\$(grep -Eo '[0-9]+(\\.[0-9]+)?' "in/yak/${sample}.yak.qv.txt" | tail -n1 || true)
  else
    QV=""
  fi
  printf "sample\\tyak_qv\\n%s\\t%s\\n" "${sample}" "\${QV:-NA}" > in/yak/yak_qv.tsv

  # Bandage image
  if [[ -f "${bandage_png}" ]]; then
    cp -f "${bandage_png}" "in/bandage/${sample}.assembly.gfa.png"
  fi

  # Write a small MultiQC config with a placeholder
  cat > multiqc_config.yaml <<'CFG'
  custom_content:
    - id: yak_qv
      section_name: "Yak QV"
      description: "Consensus quality estimate from yak (higher is better)."
      plot_type: "table"
      data_format: "tsv"
      file: "in/yak/yak_qv.tsv"
      pconfig:
        export: true

    - id: bandage_img
      section_name: "Assembly graph (Bandage)"
      description: "Rendered from the assembly GFA."
      file_format: "image"
      images:
        - "in/bandage/REPLACE_SAMPLE.assembly.gfa.png"
  CFG

  sed -i "s/REPLACE_SAMPLE/${sample}/g" multiqc_config.yaml

  multiqc -o . -c multiqc_config.yaml in
  """
}
 

