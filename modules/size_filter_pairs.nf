process SIZE_FILTER_PAIRS {

    tag { splitName }

    input:
    tuple val(splitName), file(splitSam) from resultsHicupMapper

    output:
    file("${splitName}/${splitName}_1_2.filt.sam") into resultsHicupFilter
    tuple val(splitName), file("${splitName}/*summary*.txt"), file("${splitName}/*.ditag_size_distribution") into hicupFilterReportChannel

    shell:
    '''
    mkdir !{splitName}
    filterBySize.py -i !{splitSam} \
                    --minDistance !{params.minMapDistance} \
                    -o !{splitName}/!{splitName}_1_2.filt.sam
    '''
}
