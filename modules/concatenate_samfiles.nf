process CONCATENATE_SAMFILES {

    tag { name }

    input:
    tuple val(name), file(dedupSams) from catSamInputChannel

    output:
    tuple val(name), file("${name}_1_2.dedup.sam") into resultsHicup, sam2bamChannel

    shell:
    '''
    samtools view -H !{dedupSams[0]} > !{name}_1_2.dedup.sam
    cat !{dedupSams} | grep -v '^@' >> !{name}_1_2.dedup.sam
    '''
}
