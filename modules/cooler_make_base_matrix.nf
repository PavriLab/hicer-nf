process MAKE_BASE_MATRIX {

    tag { name }

    input:
    tuple val(name), file(pairs), file(pairsIndex) from resultsPairixBaseBuilder
    file(chromSizeFile) from chromSizeChannelCooler.collect()

    output:
    tuple val(name), file("${name}/${name}_base.cool") into resultsBaseBuilder

    shell:
    '''
    mkdir -p !{name}

    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V !{chromSizeFile} > chromSizes.sort.tsv

    cooler cload pairs --assembly !{genomeName} \
                       -c1 2 -p1 3 -c2 4 -p2 5 \
                       chromSizes.sort.tsv:!{baseResolution} \
                       !{pairs} \
                       !{name}/!{name}_base.cool
    '''

}
