process HICUP_TRUNCATE_READS {

    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), file("${meta.id}_*.trunc.fastq"),  emit: reads
    tuple val(meta), file("*summary*.txt"),             emit: reports

    script:
    """
    hicup_truncater \
        --outdir ${meta.id} \
        --threads ${task.cpus} \
        --re1 ${params.re} \
        ${reads[0]} \
        ${reads[1]}
    """
}
