#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
  def s   = params.sample ?: 'Pseudo'
  def asm = params.asm     ? file(params.asm)           : null
  def lst = params.type_accessions ? file(params.type_accessions) : null
  if (!asm?.exists()) error "Assembly not found: ${params.asm}"
  if (!lst?.exists()) error "type_accessions file not found: ${params.type_accessions}"

  ch_genome = Channel.value ( file(params.genome) )

  ch_in = Channel.of( tuple(s as String, asm, lst) )
  AAI_EZAAI(ch_in, ch_genome)
}

process AAI_EZAAI {
  tag { sample }
  publishDir { "${params.outdir}/aai/${sample}" }, mode: 'copy'
  conda "${params.ezaai_env}"

  input:
    tuple val(sample), path(query_fa), path(type_list) 
    path genome_py 

  output:
    path "out/Pseudo_vs_refs.aai.tsv"
    path "out/aai.tsv"
    path "out/aai.nwk"

  shell:
  '''
  mkdir -p genomes db out

  # 1) fetch references listed (one accession per line: GCF_/GCA_)
  datasets download genome accession --inputfile "!{type_list}" \
    --include genome --filename refs.zip || true
  unzip -oq refs.zip -d refs || true
  rm -f refs.zip || true

  # 2) Stage genomes (query + labeled references) via Python helper (Py2.7 safe)
  PY="!{ params.python ?: '' }"
  if [[ -z "$PY" || ! -x "$PY" ]]; then
    PY=$(command -v python2 || command -v python || true)
  fi
  [[ -x "$PY" ]] || { echo "FATAL: no Python 2.x interpreter found"; exit 127; }

  [[ -s "!{genome_py}" ]] || { echo "FATAL: stage_genomes_py27.py not found"; exit 2; }

  "$PY" "!{genome_py}" \
    --refs-dir refs \
    --out-dir  genomes \
    --query    "!{query_fa}" \
    --sample   "!{sample}"

  # 3) EzAAI
  ezaai extract   -i genomes -o db -t !{task.cpus}
  ezaai calculate -i db -j db -o out/aai.tsv -t !{task.cpus}
  ezaai cluster   -i out/aai.tsv -o out/aai.nwk

  # 4) Pseudo-focused table for decisions
  awk 'BEGIN{FS=OFS="\\t"} NR==1 || $3 ~ /!{sample}/ || $4 ~ /!{sample}/' out/aai.tsv \
    > out/Pseudo_vs_refs.aai.tsv
  '''
}

