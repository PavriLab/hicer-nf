process MAKE_PAIRS_FILE {

    tag "$meta.id"

    publishDir path: "${params.outputDir}/${meta.id}/pairs/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "*pairs.gz*",
               saveAs: {
                   filename ->
                       if (file.endsWith(".pairs.gz")) file(file).getName()
                       else if (file.endsWith(".pairs.gz.px2")) file(file).getName()
               }

    input:
    tuple val(meta), file(alignments)
    file(genomeSizes)

    output:
    tuple val(meta), path("${meta.id}.pairs.gz*"), emit: pairs

    shell:
    '''
    sam2pairs.py !{alignments} | \
        paste - - | \
        awk 'BEGIN{ FS = "\t"; OFS = "\t" }{ print $1,$2,$3,$6,$7,$4,$8 }' > \
        !{meta.id}.pairs.tmp

    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V !{genomeSizes} > chromSizes.sort.tsv

    cooler csort \
        -c1 2 -c2 4 \
        -p1 3 -p2 5 \
        -p !{task.cpus} \
        !{meta.id}.pairs.tmp \
        chromSizes.sort.tsv

    # add generic header to make pairix compatible with juicer prefix
    echo "## pairs format v1.0" > !{meta.id}.pairs
    echo "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2" >> !{meta.id}.pairs
    zcat !{meta.id}.pairs.tmp.blksrt.gz >> !{meta.id}.pairs

    bgzip !{meta.id}.pairs
    pairix -p pairs !{meta.id}.pairs.gz
    '''
}
