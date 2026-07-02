#!/usr/bin/env nextflow

workflow {
    read_csv = Channel.fromPath("$params.index")
    | splitCsv( header: true )
    bakta(read_csv)
    
}

process bakta {
    conda '/mount/britton/Jose/conda-envs/bakta_env'
    tag "Running bakta on $sample_id (Attempt #$task.attempt)"
    publishDir "${params.outdir}/bakta/${sample_id}/", mode: 'rellink'

    input:
    // Define the input as a tuple containing:
    // - val(sample_id): a unique identifier for the file pair which acts as the basename to help separate out each file by looking at publishDir
    // - path(assembly): the path to the fasta files processed from Unicycler 
    tuple val(sample_id), path(assembly) 

    output:
    path("${sample_id}.txt")
    path("*.*")

    script:
    // Command to run Bakta with assembly fasta files processed from Unicycler 
    """
    export BAKTA_DB_PATH="/mount/britton/leap/dbs/bakta/db"
    bakta --db \$BAKTA_DB_PATH  --output ${params.outdir}/bakta/${sample_id} --force ${assembly}
    echo "Processed sample: ${sample_id}" > "${sample_id}.txt"
    echo "Genome file: ${assembly}" >> "${sample_id}.txt"
    echo "Output directory: ${params.outdir}/bakta/${sample_id}" >> "${sample_id}.txt"
    """
}


