process HICUP_MAP_READS {

    tag { splitName }
    memory = { genomeSizeType == 'large' ? 80.GB * task.attempt : 20.GB * task.attempt }

    input:
    tuple val(splitName), file(fastqTruncPairs) from resultsHicupTruncater
    file(index) from bowtie2Index.collect()

    output:
    tuple val(splitName), file("${splitName}/${splitName}_1_2.pair.sam") into resultsHicupMapper
    tuple val(splitName), file("${splitName}/*summary*.txt") into hicupMapperReportChannel

    shell:
    '''
    mkdir !{splitName}
    hicup_mapper --outdir !{splitName} \
                 --threads !{task.cpus} \
                 --format Sanger \
                 --index !{index}/!{bwt2_base} \
                 --bowtie2 $(which bowtie2) \
                 !{fastqTruncPairs[0]} \
                 !{fastqTruncPairs[1]}
    '''
}
