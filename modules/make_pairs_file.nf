process MAKE_PAIRS_FILE {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/pairs/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*pairs.gz*",
               saveAs: { filename ->
                             if (filename.endsWith(".pairs.gz")) file(filename).getName()
                             else if (filename.endsWith(".pairs.gz.px2")) file(filename).getName()
                       }

    input:
    tuple val(name), file(sam) from resultsHicup
    file(chromSizeFile) from chromSizeChannelPairix.collect()

    output:
    tuple val(name), file("${name}/${name}.pairs.gz"), file("${name}/${name}.pairs.gz.px2") into resultsPairixBaseBuilder, resultsPairixJuicer

    shell:
    '''
    mkdir -p !{name}

    sam2pairs.py !{sam} | \
        paste - - | \
        awk 'BEGIN{ FS = "\t"; OFS = "\t" }{ print $1,$2,$3,$6,$7,$4,$8 }' > \
        !{name}/!{name}.pairs.tmp

    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V !{chromSizeFile} > chromSizes.sort.tsv

    cooler csort -c1 2 -c2 4 \
                 -p1 3 -p2 5 \
                 -p !{task.cpus} \
                 !{name}/!{name}.pairs.tmp \
                 chromSizes.sort.tsv

    # add generic header to make pairix compatible with juicer prefix
    echo "## pairs format v1.0" > !{name}/!{name}.pairs
    echo "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2" >> !{name}/!{name}.pairs
    zcat !{name}/!{name}.pairs.tmp.blksrt.gz >> !{name}/!{name}.pairs

    bgzip !{name}/!{name}.pairs
    pairix -p pairs !{name}/!{name}.pairs.gz
    '''
}
