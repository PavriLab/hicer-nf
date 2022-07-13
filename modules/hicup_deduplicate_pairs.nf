process HICUP_DEDUPLICATE_PAIRS {

    tag { resplitName }

    input:
    file sam from resultsResplit.flatten()

    output:
    file("${resplitName}/${resplitName}_1_2.dedup.sam") into resultsHicupDeduplicate
    tuple val(resplitName), file("${resplitName}/*summary*.txt") into hicupDeduplicatorReportChannel

    shell:
    resplitName = (sam.getName() - ~/(_1_2\.filt\.sam)?$/)

    '''
    mkdir !{resplitName}
    hicup_deduplicator --outdir !{resplitName} \
                       !{sam}

    if [ ! !{params.re} ]
    then
        getCisFractions.py -i !{resplitName}/!{resplitName}_1_2.dedup.sam \
                           -r !{resplitName}/*summary*.txt
    fi
    '''
}
