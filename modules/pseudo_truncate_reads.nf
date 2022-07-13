process PSEUDO_TRUNCATE_READS {

    tag { splitName }

    input:
    tuple val(splitName), file(fastqSplitPairs) from truncaterInputChannel

    output:
    tuple val(splitName), file("${splitName}/${splitName}_*.trunc.fastq") into resultsHicupTruncater
    tuple val(splitName), file("${splitName}/*summary*.txt") into hicupTruncaterReportChannel

    shell:
    '''
    mkdir !{splitName}
    dummyReportGenerator.py -1 !{fastqSplitPairs[0]} \
                            -2 !{fastqSplitPairs[1]} \
                            -o !{splitName}/!{splitName}.dummy_truncater_summary.txt

    cp !{fastqSplitPairs[0]} !{splitName}/!{splitName}_1.trunc.fastq
    cp !{fastqSplitPairs[1]} !{splitName}/!{splitName}_2.trunc.fastq
    '''
}
