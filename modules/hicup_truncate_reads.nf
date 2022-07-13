process HICUP_TRUNCATE_READS {

    tag { splitName }

    input:
    tuple val(splitName), file(fastqSplitPairs) from truncaterInputChannel

    output:
    tuple val(splitName), file("${splitName}/${splitName}_*.trunc.fastq") into resultsHicupTruncater
    tuple val(splitName), file("${splitName}/*summary*.txt") into hicupTruncaterReportChannel

    shell:
    '''
    mkdir !{splitName}
    hicup_truncater --outdir !{splitName} \
                    --threads !{task.cpus} \
                    --re1 !{params.re} \
                    !{fastqSplitPairs[0]} \
                    !{fastqSplitPairs[1]}
    '''
}
