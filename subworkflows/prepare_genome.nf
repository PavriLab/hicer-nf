include { GUNZIP as GUNZIP_FASTA } from '../modules/gunzip.nf'
include { XML_TO_TSV             } from '../modules/xml_to_tsv.nf'
include { DIGEST_GENOME          } from '../modules/hicup_digest_genome.nf'
include { BOWTIE2_BUILD_INDEX    } from '../modules/bowtie2_build_index.nf'

workflow PREPARE_GENOME {
    take:
    prepare_genome_for_tools
    dynamic_params

    main:

    // Uncompress genome fasta file if required
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( dynamic_params.genomeFasta )

    } else {
        ch_fasta = file( dynamic_params.genomeFasta )
    }

    if ("bowtie2" in prepare_genome_for_tools) {
        ch_bowtie2_index = BOWTIE2_BUILD_INDEX (
            ch_fasta,
            dynamic_params.genomeSizeType
        )

    } else {
        ch_bowtie2_index = file( dynamic_params.bowtie2Index )
    }

    ch_genome_digest = Channel.empty()
    if ("hicup" in prepare_genome_for_tools) {
        ch_genome_digest = DIGEST_GENOME (
            ch_fasta,
            dynamic_params.re
        )
    }

    if ("chromSizes" in prepare_genome_for_tools) {
        ch_genome_sizes = XML_TO_TSV ( dynamic_params.genomeSizes )

    } else {
        ch_genome_sizes = file( dynamic_params.genomeSizes )
    }

    emit:
    index      = ch_bowtie2_index
    digest     = ch_genome_digest
    sizes      = ch_genome_sizes

}
