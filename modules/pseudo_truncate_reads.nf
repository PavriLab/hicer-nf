process PSEUDO_TRUNCATE_READS {

    tag "$meta.id"

    input:
    tuple val(meta), file(reads)

    output:
    tuple val(meta), file("${meta.id}_*.trunc.fastq"),  emit: reads
    tuple val(meta), file("*summary*.txt"),             emit: reports

    script:
    """
    dummyReportGenerator.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o ${meta.id}.dummy_truncater_summary.txt

    cp ${reads[0]} ${meta.id}_1.trunc.fastq
    cp ${reads[1]} ${meta.id}_2.trunc.fastq
    """
}
