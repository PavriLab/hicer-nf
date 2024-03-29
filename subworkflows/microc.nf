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
    bowtie2Index
    genomeSizeType

    main:
    PSEUDO_TRUNCATE_READS ( ch_fastq )

    HICUP_MAP_READS (
        PSEUDO_TRUNCATE_READS.out.reads,
        bowtie2Index,
        genomeSizeType
    )

    SIZE_FILTER_PAIRS ( HICUP_MAP_READS.out.alignments )
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
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = file.name.toString() - ~/(_1_2\.filt\.sam)?$/
                [ meta_clone, file ]
        }
        .set { ch_resplit_pairs }

    HICUP_DEDUPLICATE_PAIRS ( ch_resplit_pairs )
        .alignments
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = file.name.toString() - ~/(_[a-zA-Z0-9]+_1_2\.dedup\.sam)?$/
                [ meta_clone, file ]
        }
        .groupTuple(by: [0])
        .set { ch_cat_sam }

    CAT_SAM ( ch_cat_sam )

    GENERATE_HICUP_REPORT (
        PSEUDO_TRUNCATE_READS.out.reports,
        HICUP_MAP_READS.out.reports,
        SIZE_FILTER_PAIRS.out.reports,
        HICUP_DEDUPLICATE_PAIRS.out.reports
    )

    emit:
    alignments  = CAT_SAM.out.alignments
    multiqc     = GENERATE_HICUP_REPORT.out.multiqc
}
