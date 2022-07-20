include { HICUP_TRUNCATE_READS    } from '../modules/hicup_truncate_reads.nf'
include { HICUP_MAP_READS         } from '../modules/hicup_map_reads.nf'
include { HICUP_FILTER_PAIRS      } from '../modules/hicup_filter_pairs.nf'
include { RESPLIT_FILTERED_PAIRS  } from '../modules/resplit_filtered_pairs.nf'
include { HICUP_DEDUPLICATE_PAIRS } from '../modules/hicup_deduplicate_pairs.nf'
include { CAT_SAM                 } from '../modules/cat_sam.nf'
include { GENERATE_HICUP_REPORT   } from '../subworkflows/generate_hicup_report.nf'

workflow HIC {
    take:
    ch_fastq
    bowtie2Index
    genomeDigest
    genomeSizeType

    main:
    HICUP_TRUNCATE_READS ( ch_fastq )

    HICUP_MAP_READS (
        HICUP_TRUNCATE_READS.out.reads,
        bowtie2Index,
        genomeSizeType
    )

    HICUP_FILTER_PAIRS (
        HICUP_MAP_READS.out.alignments,
        genomeDigest
    )

    HICUP_FILTER_PAIRS
        .alignments
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = file.name.toString() - ~/(_[a-z]{4}_1_2\.filt\.sam)?$/
                [ meta_clone, file ]
        }
        .groupTuple (by: [0])
        .set { ch_filtered_pairs }

    RESPLIT_FILTERED_PAIRS ( ch_filtered_pairs )
        .alignments
        .map { WorkflowHicer.distributeMetaSingle( it ) }
        .flatten ()
        .collate ( 2 )
        .set { ch_resplit_pairs }

    HICUP_DEDUPLICATE_PAIRS ( ch_resplit_pairs )
        .alignments
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
