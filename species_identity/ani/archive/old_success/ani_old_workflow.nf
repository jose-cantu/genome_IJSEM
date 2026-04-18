#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------- ANI-only entry (uses existing BUILD_FASTANI_PANEL & FASTANI_BATCH) ----------
workflow {

  // required params: --sample, --asm, --class_summary
  def s        = params.sample ?: null
  def asmFile  = params.asm ? file(params.asm) : null
  def sumFile  = params.class_summary ? file(params.class_summary) : null

  if (!s)                          error "Provide --sample"
  if (!asmFile?.exists())          error "Assembly not found: ${params.asm}"
  if (!sumFile?.exists())          error "GTDB summary not found: ${params.class_summary}"

  // Channels
  asm_ch   = Channel.of( tuple(s as String, asmFile) )
  class_ch = Channel.of( tuple(s as String, sumFile) )

  // Build refs.txt from GTDB summary, then run FastANI (process code unchanged)
  panel_ch   = BUILD_FASTANI_PANEL(class_ch)
  fastani_in = asm_ch.join(panel_ch).map { sample, asm_fa, refs, meta_dir -> tuple(sample, asm_fa, refs) }
  fastani_out = FASTANI_BATCH(fastani_in)
  
  // aannotatie with TYPE flags (joining fastANI + META)
  annotate_in = fastani_out.join(panel_ch).map { smp, fastani_tsv, smp2, refs, meta_dir -> tuple(smp, fastani_tsv, meta_dir) }
  FASTANI_ANNOTATE_TYPES(annotate_in) 
}

// 12a) Build FastANI reference panel from GTDB summary
process BUILD_FASTANI_PANEL {
  tag { sample }
  publishDir { "${params.outdir}/fastani_panel/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/ncbi_datasets"

  input:
    tuple val(sample), path(class_summary)

  output:
    tuple val(sample), path("refs.txt"), path("meta")

  // NOTE: use triple double quotes + Groovy defaults inside !{ }
  shell:
  '''
  
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

// Post-process: join fastANI with NCBI metadata -> fastani_annotated.tsv
process FASTANI_ANNOTATE_TYPES {
  tag { sample }
  // publish alongside fastani.tsv
  publishDir { "${params.outdir}/fastani/${sample}" }, mode: 'copy'
  // needs `datasets`, `jq` (same env you used for panel)
  conda "/mount/britton/Jose/conda-envs/ncbi_datasets"

  input:
  // from FASTANI_BATCH + BUILD_FASTANI_PANEL
  tuple val(sample), path(fastani_tsv), path(ref_list)

  output:
  path("fastani_annotated.tsv")

  shell:
  ''' 
  # 1) Build meta.tsv from JSONL
  jq -r '[
      (.accession // "NA"),
      (.organism.organism_name // "NA"),
      (.assembly_info.relation_to_type_material // "NA"),
      (.assembly_info.assembly_level // "NA"),
      (.assembly_info.assembly_status // "NA")
    ] | @tsv' "\${meta_dir}"/*.jsonl 2>/dev/null | sort -k1,1 > meta.tsv || : > meta.tsv

  # 2) Reduce FastANI to accession + ANI + AF  (AF = aligned/total fragments)
  awk -F '\\t' '{
     n=split(\$2,a,"/"); acc=a[n]; sub(/\\.fna$/,"",acc);
     af=(NF>=5 && \$5>0) ? \$4/\$5 : "";
     printf "%s\\t%.4f\\t%s\\n", acc, \$3, (af==""? "": sprintf("%.3f",af));
  }' "\${ani_tsv}" | sort -k1,1 > ani.tsv

  # 3) Join -> annotated + TYPE flag
  join -t $'\\t' -1 1 -2 1 ani.tsv meta.tsv \
  | awk -F'\\t' 'BEGIN{
      OFS="\\t";
      print "accession","ANI","AF","organism","relation_to_type_material","assembly_level","status","TYPE_flag"
    }{
      t=tolower(\$5)~/(^| )type( |$)/ ? "TYPE" : "non-type";
      print \$1,\$2,\$3,\$4,\$5,\$6,\$7,t
    }' > fastani_annotated.tsv
  '''
}

