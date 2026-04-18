#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------- ANI-only entry (uses BUILD_FASTANI_PANEL & FASTANI_BATCH) ----------
workflow {

  // required params: --sample, --asm, --class_summary
  def sample   = params.sample ?: null
  def asmFile  = params.asm ? file(params.asm) : null
  def sumFile  = params.class_summary ? file(params.class_summary) : null

  if (!s)                 error "Provide --sample"
  if (!asmFile?.exists()) error "Assembly not found: ${params.asm}"
  if (!sumFile?.exists()) error "GTDB summary not found: ${params.class_summary}"

  // Channels
  asm_ch   = Channel.of( tuple(s as String, asmFile) )
  class_ch = Channel.of( tuple(s as String, sumFile) )
  CH_TOK = Channel.value ( file(params.gtdb_tokens) ) 

  // Build refs.txt (NCBI type-first; fall back to GTDB neighbours), return meta/ for annotation
  panel_ch   = BUILD_FASTANI_PANEL(class_ch, CH_TOK)

  // Prepare FastANI input
  fastani_in  = asm_ch.join(panel_ch).map { sample, asm_fa, refs, meta_dir -> tuple(sample, asm_fa, refs) }
  fastani_out = FASTANI_BATCH(fastani_in)

  // Annotate with TYPE flags (join FastANI + NCBI meta/)
  annotate_in = fastani_out.join(panel_ch).map { smp, fastani_tsv, smp2, refs, meta_dir -> tuple(smp, fastani_tsv, meta_dir) }
  FASTANI_ANNOTATE_TYPES(annotate_in)
}

// ---------- Build FastANI reference panel (TYPE-first) ----------
process BUILD_FASTANI_PANEL {
  tag { sample }
  publishDir { "${params.outdir}/fastani_panel/${sample}" }, mode: 'copy'
  conda "/mount/britton/Jose/conda-envs/ncbi_datasets"

  input:
    tuple val(sample), path(class_summary)
    path TOK_PY 

  // Return both refs.txt and the metadata directory for annotation
  output:
    tuple val(sample), path("refs.txt"), path("meta")

  shell:
  '''
  # 0) Prepare output directories
  mkdir -p panel meta

  # 1) Try NCBI “type material” discovery by taxon first
  if [[ "!{ params.fastani_ncbi_from_type ?: true }" == "true" ]]; then
    : > accessions.type.txt
    # Extract genus (g__ token) and family from the GTDB classification column
    PY=""
    for cand in "!{params.python ?: ''}" python3 python /usr/bin/python3 /usr/bin/python /usr/local/bin/python3 /usr/local/bin/python; do
     [[ -z "$cand" ]] && continue
     if command -v "$cand" >/dev/null 2>&1; then PY="$(command -v "$cand")"; break; fi
     [[ -x "$cand" ]] && { PY="$cand"; break; }
    done
    [[ -z "$PY" ]] && { echo "FATAL: no Python interpreter on PATH (tried params.python, python3/python)"; exit 127; }
    echo "USING PY: $PY" >&2
    read iso_genus iso_family < <( "$PY" "!{params.gtdb_tokens}" "!{class_summary}" ) # parse tokens once ......  

    echo "TYPE-SCAN genus=${iso_genus:-NA} family=${iso_family:-NA}" >&2

    if [[ -n "${iso_genus:-}" ]]; then
      # Query NCBI for type-material assemblies in this genus (limit by fastani_max_refs)
      datasets summary genome taxon "${iso_genus}" --from-type --report genome --as-json-lines | jq -r '.. | objects | .accession? // .assembly?.accession? // empty' | grep -E '^GC[AF]_' | sort -u | head -n !{ params.fastani_max_refs ?: 20 } > accessions.type.txt || true
    fi
     
    # If genus yielded none, try family-level type material (e.g., Oscillospiraceae)   
    if [[ ! -s accessions.type.txt && -n "${iso_family:-}" ]]; then
        datasets summary genome taxon "${iso_family}" --from-type --report genome --as-json-lines | jq -r '.. | objects | .accession? // .assembly?.accession? // empty' | grep -E '^GC[AF]_' | sort -u | head -n !{ params.fastani_max_refs ?: 20 } > accessions.type.txt || true 
    fi 

    echo "TYPE-SCAN found $( [ -s accessions.type.txt ] && wc -l < accessions.type.txt || echo 0 ) candidates" >&2
    # If we found any type-material accessions, download and use them
    if [[ -s accessions.type.txt ]]; then
      while read -r acc; do
        [[ -z "$acc" ]] && continue
        datasets summary genome accession "$acc" --as-json-lines > "meta/${acc}.jsonl" || true
        datasets download genome accession "$acc" --include genome --filename "${acc}.zip" || continue
        unzip -oq "${acc}.zip" -d "${acc}" || true
        f=$(find "${acc}/ncbi_dataset/data" -type f \\( -name "*_genomic.fna" -o -name "*.fna" \\) | head -n1 || true)
        [[ -n "$f" ]] && cp -f "$f" "panel/${acc}.fna" 
        rm -rf "${acc}.zip" "${acc}"
      done < accessions.type.txt
      if compgen -G "panel/*.fna" > /dev/null; then
        realpath panel/*.fna > refs.txt
        exit 0
      fi
    fi
  fi

  # 2) Fall back to GTDB neighbours (original logic)
  awk -F '\\t' '
    function emit_token(tok,    p,acc) {
      p = index(tok,"GCF_"); if (!p) p = index(tok,"GCA_");
      if (p) {
        acc = substr(tok,p)
        sp = index(acc," "); cm = index(acc,","); sc = index(acc,";"); cut = 0
        if (sp && (!cut || sp<cut)) cut = sp
        if (cm && (!cut || cm<cut)) cut = cm
        if (sc && (!cut || sc<cut)) cut = sc
        if (cut) acc = substr(acc,1,cut-1)
        if (index(acc,".")>0) print acc
      }
    }
    NR==1{
      for (i=1;i<=NF;i++) h[$i]=i
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
      if (pr  && $pr!=""  && $pr!="NA")  emit_token($pr)
      if (orr && $orr!="" && $orr!="NA") {
          s = $orr; gsub(",", " ", s); gsub(";", " ", s)
          n = split(s,a," ")
          for (i=1;i<=n;i++) if (a[i]!="") emit_token(a[i])
      }
    }
  ' "!{class_summary}" | sort -u | head -n !{ params.fastani_max_refs ?: 20 } > accessions.txt

  # Download the GTDB neighbour accessions
  while read -r acc; do
    datasets summary genome accession "$acc" --as-json-lines > "meta/${acc}.jsonl" || true
    datasets download genome accession "$acc" --include genome --filename "${acc}.zip" || continue
    unzip -oq "${acc}.zip" -d "${acc}" || true
    f=$(find "${acc}/ncbi_dataset/data" -type f \\( -name "*_genomic.fna" -o -name "*.fna" \\) | head -n1 || true)
    if [[ -n "$f" ]]; then
      cp -f "$f" "panel/${acc}.fna"
    fi
    rm -rf "${acc}.zip" "${acc}"
  done < accessions.txt

  # 3) Build TYPE and ALL lists, then enforce require/prefer policies
  : > refs.all.txt
  : > refs.type.txt
  if compgen -G "panel/*.fna" > /dev/null; then
    while read -r acc; do
      f="panel/${acc}.fna"; [[ -s "$f" ]] || continue
      p=$(realpath "$f"); echo "$p" >> refs.all.txt
      if jq -er '.assembly_info.relation_to_type_material // empty | ascii_downcase | test("type")' "meta/${acc}.jsonl" > /dev/null 2>&1; then
        echo "$p" >> refs.type.txt
      fi
    done < accessions.txt
  fi

  if [[ "!{ params.fastani_require_type ?: false }" == "true" ]]; then
    if [[ -s refs.type.txt ]]; then
      sort -u refs.type.txt > refs.txt
    else
      echo "ERROR: No 'type material' genomes found among GTDB neighbours." >&2
      : > refs.txt
      exit 1
    fi
  elif [[ "!{ params.fastani_prefer_type ?: false }" == "true" && -s refs.type.txt ]]; then
    sort -u refs.type.txt > refs.txt
  elif [[ -s refs.all.txt ]]; then
    sort -u refs.all.txt > refs.txt
  else
    : > refs.txt
  fi
  '''
}

// ---------- FastANI ----------
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

// ---------- Annotate FastANI with NCBI metadata ----------
process FASTANI_ANNOTATE_TYPES {
  tag { sample }
  publishDir { "${params.outdir}/fastani/${sample}" }, mode: 'copy'
  conda "/home/britton/Jose/conda-envs/fastani"

  input:
    tuple val(sample), path(fastani_tsv), path(meta_dir)

  output:
    path("fastani_annotated.tsv")

  shell:
  '''
  set -euo pipefail
  shopt -s nullglob

  # 1) Build meta.tsv from JSONL (accession, organism, type relation, assembly level, status, length, GC%)
  jq -r '[
      (.accession // "NA"),
      (.organism.organism_name // "NA"),
      (.assembly_info.relation_to_type_material // "NA"),
      (.assembly_info.assembly_level // "NA"),
      (.assembly_info.assembly_status // "NA"),
      (.assembly_stats.total_sequence_length // "NA"),
      (.assembly_stats.gc_percent // "NA")
    ] | @tsv' "!{meta_dir}"/*.jsonl 2>/dev/null | sort -k1,1 > meta.tsv || : > meta.tsv

  # 2) Reduce FastANI to accession + ANI + AF  (AF = aligned/total fragments)
  awk -F '\\t' 'NR>1{
     n=split($2,a,"/"); acc=a[n]; sub(/\\.fna$/,"",acc);
     af=(NF>=5 && $5>0) ? $4/$5 : "";
     printf "%s\\t%.4f\\t%s\\n", acc, $3, (af==""? "": sprintf("%.3f",af));
  }' "!{fastani_tsv}" | sort -k1,1 > ani.tsv

  # 3) Join -> annotated + TYPE flag (keep ANI rows even if metadata missing)
  join -t $'\\t' -a 1 -e 'NA' -o auto -1 1 -2 1 ani.tsv meta.tsv \
  | awk -F'\\t' 'BEGIN{
      OFS="\\t";
      print "accession","ANI","AF","organism","relation_to_type_material",
            "assembly_level","status","genome_bp","gc_percent","TYPE_flag"
    }{
      t=tolower($5)~/(^| )type( |$)/ ? "TYPE" : "non-type";
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,t
    }' > fastani_annotated.tsv
  '''
}

