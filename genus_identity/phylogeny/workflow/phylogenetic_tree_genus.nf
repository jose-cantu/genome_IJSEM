#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {

  def csvPath = params.index?.trim()
  def csvFile = file(csvPath)
  if( !csvFile.exists() )
      error "CSV not found on this node: ${csvFile}"

  // Create a single-value channel carrying the CSV file, then split rows
  reads_ch = Channel.value(csvFile)
      .splitCsv(header: true)
      // drop bad/incomplete lines early
      .filter { row -> row.sample && row.ont && (!params.sample || row.sample == params.sample) }
      .map { row -> tuple(row.sample as String, file(row.ont)) }

  // QC and Trimming 
  ont_ch = PORECHOP_ABI(reads_ch)
  
  // Flye assembly  
  asm_ch   = FLYE(ont_ch)

  // validate params.type_accession before download 
  def accFile = file(params.type_accession ?: "")
  if ( !accFile.exists() )
       error "type_accession file not found: ${params.type_accession}" 
 
  // download curated accession list from TYGS and LPSN API 
  accession_ch = Channel.value(accFile)
  refs_ch = DOWNLOAD_REFS(accession_ch) 

  // QC targets: Pseudo assembly + all downloaded reference genomes
  isolate_qc_ch = asm_ch.map { sample, asm_fa ->
      tuple("ISO_${sample}" as String, asm_fa)
  }

  ref_qc_ch = refs_ch.flatMap { refs_dir ->
      refs_dir.toFile().listFiles()
          .findAll { it.name.endsWith('.fna') }
          .collect { f ->
              def id = f.name.replaceFirst(/\.fna$/, '')
              tuple(id as String, file(f.toString()))
          }
  }

  qc_targets_ch = isolate_qc_ch.mix(ref_qc_ch)

  quast_qc   = QUAST_QC(qc_targets_ch)
  checkm1_qc = CHECKM1_QC(qc_targets_ch)

  panel_ch = asm_ch.combine(refs_ch)  

  denovo_ch = GTDB_DENOVO_ALIGN(panel_ch) 
  IQTREE2(denovo_ch) 

}


/* --------------- PORECHOP_ABI Process ---------------- */

// ONT adapter trim (ab initio)
process PORECHOP_ABI {
  conda '/mount/britton/Jose/conda-envs/porechop_abi'
  tag { sample }
  publishDir { "${params.outdir}/porechop_abi/${sample}" }, mode: 'copy'

  input:
  tuple val(sample), path(ont)

  output:
  tuple val(sample),
        path("${sample}.ont.trim.fastq.gz")

  script:
  """
  # Porechop_ABI: infer ONT adapters directly from reads, then trim
  porechop_abi -abi \
    -i "${ont}" \
    -o "${sample}.ont.trim.fastq.gz" \
    -t ${task.cpus}
  """
}


// Flye (long read first assembly approach dnaapler rotate replicons to dnaA/repA at next step)
process FLYE {
  tag { sample }
  publishDir { "${params.outdir}/Dnaapler/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/mamba/envs/autocycler_env" 

  input: 
  tuple val(sample), path(ont_trim) 
 
  output: 
  tuple val(sample), path("assembly.fasta") 

  script:
  """
  flye --nano-raw "${ont_trim}" --out-dir flye --genome-size 2.6m --asm-coverage 100 --threads ${task.cpus}
  cp -f flye/assembly.fasta assembly.fasta 
  """ 
}

process QUAST_QC {
  tag { id }
  publishDir { "${params.outdir}/qc/quast/${id}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/quast"

  input:
  tuple val(id), path(fa)

  output:
  tuple val(id), path("quast")

  script:
  """
  quast.py "${fa}" -o quast -t ${task.cpus}
  """
}

process CHECKM1_QC {
  tag { id }
  publishDir { "${params.outdir}/qc/checkm1/${id}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/checkm"

  input:
  tuple val(id), path(fa)

  output:
  tuple val(id), path("checkm1"), path("checkm.qa.tsv")

  script:
  """
  mkdir -p bins
  cp -f "${fa}" "bins/${id}.fna"

  checkm lineage_wf -x fna -t ${task.cpus} bins checkm1
  checkm qa -o 2 checkm1/lineage.ms checkm1 > checkm.qa.tsv
  """
}

process DOWNLOAD_REFS {
  tag "refs" 
  publishDir { "${params.outdir}/refs_genomes" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/ncbi_datasets"

  input:
  path(accession_file) 

  output:
  path("refs_genomes")

  script:
  """
  mkdir -p refs_genomes 
  datasets download genome accession --inputfile "${accession_file}" --include genome --filename refs.zip 
  unzip -q refs.zip -d refs_ncbi

  # Normalize filenames: Extract <accession> from <accession>_<assembly>_genomic.fna
  find refs_ncbi -type f -name "*_genomic.fna" -print0 |
    while IFS= read -r -d '' f; do
      bn=\$(basename "\$f")        # e.g., GCA_002834225.1_ASM283422v1_genomic.fna
      acc=\$(printf "%s" "\$bn" | cut -d'_' -f1-2)  # GCA_002834225.1
      cp -f "\$f" "refs_genomes/\${acc}.fna"
    done
  
  echo "Staged genomes (head):" >&2
  ls -1 refs_genomes | head -n 10 >&2

  n=\$(ls -1 refs_genomes/*.fna 2>/dev/null | wc -l)
  echo "Total staged refs: \$n" >&2

  # Quick sanity check to make sure I downloaded it correctly 
  echo "Downloaded genomes to check refs" >&2  
  n=\$(ls -1 refs_genomes/*.fna 2>/dev/null | wc -l)
  echo "Downloaded refs \$n" >&2
  [[ "\$n" -ge 1 ]] || { echo "FATAL no reference genomes downloaded sooo stoping here and fix it"; exit 2; } 
  """
}

// gtdb devo alignment tree build de-novo marker MSA 
process GTDB_DENOVO_ALIGN { 
  tag { sample } 
  publishDir { "${params.outdir}/trees/gtdb_denovo/${sample}" }, mode: 'copy' 
  conda "/mount/britton/Jose/conda-envs/gtdbtk" 
  
  input: 
  tuple val(sample), path(asm_fa), path(refs_dir) 
  
  output: 
  tuple val(sample),
  path("user_msa.fasta.gz"), // masked concatenated MSA 
  path("fasttree.tree") // quick tree reference 
  
  script:
  """ 
  export GTDBTK_DATA_PATH="${params.gtdbtk_db}" 
   # Build a renamed genome panel so GTDB-Tk IDs never collide with GTDB reps
  rm -rf genomes_renamed
  mkdir -p genomes_renamed

  # 1) Copy isolate (give it a non-accession ID too, for consistency)
  iso_id="ISO_${sample}"
  cp -f "${asm_fa}" "genomes_renamed/\${iso_id}.fna"

  # 2) Create taxonomy + optional id map
  : > custom_taxonomy.tsv
  : > id_map.tsv

  # taxonomy line for isolate (ingroup)
  printf "%s\td__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__;s__\n" \
    "\${iso_id}" >> custom_taxonomy.tsv

  # 3) Add curated refs, rename them, and label outgroup in taxonomy
  out_hits=0

  for f in "${refs_dir}"/*.fna; do
    old_id=\$(basename "\$f" .fna)          # e.g. GCA_002834225.1
    new_id="REF_\${old_id}"               # avoid GTDB collision

    cp -f "\$f" "genomes_renamed/\${new_id}.fna"
    printf "%s\t%s\n" "\${new_id}" "\${old_id}" >> id_map.tsv

    if [[ "\${old_id}" == "${params.outgroup_acc}"* ]]; then
      printf "%s\td__Bacteria;p__OUTGROUP;c__;o__;f__;g__;s__\n" "\${new_id}" >> custom_taxonomy.tsv
      out_hits=\$((out_hits+1))
    else
      printf "%s\td__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__;s__\n" \
        "\${new_id}" >> custom_taxonomy.tsv
    fi
  done

  # hard fail if outgroup not found (or found multiple times)
  [[ "\${out_hits}" -eq 1 ]] || {
    echo "FATAL: expected exactly 1 outgroup matching '${params.outgroup_acc}*' but found \${out_hits}" >&2
    echo "DEBUG: available ref IDs:" >&2
    ls -1 "${refs_dir}"/*.fna | sed 's#.*/##; s/\\.fna\$//' >&2
    exit 2
  }

  # Here all I need is type strains + type species anchors + isolate + outgroup and avoid GTDBTK adding references to explode taxon count and make IQ-TREE slow this is more than enough I filtered to family group in reference only....
  gtdbtk de_novo_wf --genome_dir genomes_renamed --bacteria --skip_gtdb_refs --custom_taxonomy_file custom_taxonomy.tsv --outgroup_taxon p__OUTGROUP --out_dir denovo --cpus ${task.cpus} -x fna 
  
  M=\$(find denovo -type f -name "*msa*.fasta.gz" -print -quit || true) 
  T=\$(find denovo -type f -name "*.tree" -print -quit || true) 
  [[ -n "\$M" ]] && cp -f "\$M" user_msa.fasta.gz || { echo "MSA not found"; exit 2; } 
  [[ -n "\$T" ]] && cp -f "\$T" fasttree.tree || printf "(%s);\n" "${sample}" > fasttree.tree 
  
  """ 
    
} 

// reinfer tree w IQTREE2 (decided not to use ModelFinder to computational expensive used a base model already known for  + UFBoost +SH-aLRT) 
process IQTREE2 { 
  tag { sample } 
  publishDir { "${params.outdir}/trees/iqtree/${sample}" }, mode: 'copy' 
  conda "/mount/britton/Jose/conda-envs/iqtree" 
  
  input:
  tuple val(sample), path(msa_gz), path(ftree) 
  
  output:
  tuple val(sample),
  path("iqtree/${sample}.treefile"),
  path("iqtree/${sample}.iqtree"),
  path("iqtree/${sample}.log") 
  
  script:
  """
  mkdir -p iqtree 
  gunzip -c "${msa_gz}" > msa.fasta
  
  ntax=\$(grep -c '^>' msa.fasta || true) 
  echo "Ntax in MSA is: \$ntax" >&2 
  [[ "\$ntax" -ge 4 ]] || { echo "FATAL you need here >=4 taxa José for meaningful support go check your logic and fix this issue"; exit 2; }
  
  iqtree2 -s msa.fasta -m MFP -B 1000 --alrt 1000 -T ${task.cpus} --seed 42 -pre "iqtree/${sample}" 
  """ 
}



