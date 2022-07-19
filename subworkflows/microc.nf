include { PSEUDO_TRUNCATE_READS   } from '../modules/pseudo_truncate_reads.nf'
include { HICUP_MAP_READS         } from '../modules/hicup_map_reads.nf'
include { SIZE_FILTER_PAIRS       } from '../modules/size_filter_pairs.nf'
include { RESPLIT_FILTERED_PAIRS  } from '../modules/resplit_filtered_pairs.nf'
include { HICUP_DEDUPLICATE_PAIRS } from '../modules/hicup_deduplicate_pairs.nf'
include { CAT_SAM                 } from '../modules/cat_sam.nf'
include { GENERATE_HICUP_REPORT   } from '../subworkflows/generate_hicup_report.nf'

workflow MICROC {
    take:
    ch_fastq
    ch_genome
    genomeSizeType

    main:
    PSEUDO_TRUNCATE_READS ( ch_fastq )

    HICUP_MAP_READS (
        HICUP_TRUNCATE_READS.out.reads,
        ch_genome.index,
        genomeSizeType
    )

    SIZE_FILTER_PAIRS (
        HICUP_MAP_READS.out.alignments,
        ch_genome.digest
    )

    HICUP_FILTER_PAIRS.out.alignments
        .map {
        meta, file ->
              def meta_clone = meta.clone()
              meta_clone.id = file.name.toString() - ~/(_[a-z]{4}_1_2\.filt\.sam)?$/
              [ meta, file ]
        }
        .groupTuple ()
        .set { ch_filtered_pairs }

    RESPLIT_FILTERED_PAIRS ( ch_filtered_pairs )

    RESPLIT_FILTERED_PAIRS.out.alignments
        .map { WorkflowHicer.distributeMeta( it ) }
        .flatten()
        .set { ch_resplit_pairs }

    HICUP_DEDUPLICATE_PAIRS ( ch_resplit_pairs )

        HICUP_DEDUPLICATE_PAIRS.out.alignments
        .groupTuple(by: [0])
        .set { ch_cat_sam }

    CAT_SAM ( ch_cat_sam )

    GENERATE_HICUP_REPORT (
        HICUP_TRUNCATE_READS.out.reports,
        HICUP_MAP_READS.out.reports,
        HICUP_FILTER_PAIRS.out.reports,
        HICUP_DEDUPLICATE_PAIRS.out.reports
    )

    emit:
    alignments  = CAT_SAM.out.alignments
    multiqc     = GENERATE_HICUP_REPORT.out.multiqc
}
