// this subworkflow was adapted from nf-core/rnaseq
workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    Channel
        .fromPath ( samplesheet )
        .splitCsv ( header:true, sep:'\t' )
        .view()
        // .map { create_fastq_channel(it) }
        // .set { reads }

    emit:
    reads                                     // channel: [ val(sample), [ reads ] ]
}

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample

    if(!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if(!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_2}"
    }

    return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
}
