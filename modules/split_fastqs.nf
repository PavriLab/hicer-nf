process SPLIT_FASTQS {

    tag { name }

    input:
    tuple val(name), file(fastq1), file(fastq2) from resultsTrimming

    output:
    file("${name}/${name}_*") into resultsSplitting

    shell:
    numberOfLinesPerSplit = params.readsPerSplit.toInteger() * 4
    '''
    mkdir !{name}

    # ampersand at the end of the line lets linux execute the commands in parallel
    zcat !{fastq1} | \
    split -l !{numberOfLinesPerSplit} -a 4 --additional-suffix _1.fq - !{name}/!{name}_ &

    zcat !{fastq2} | \
    split -l !{numberOfLinesPerSplit} -a 4 --additional-suffix _2.fq - !{name}/!{name}_
    '''
}
