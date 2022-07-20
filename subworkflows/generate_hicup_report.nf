include { HICUP_GENERATE_REPORT } from '../modules/hicup_generate_report.nf'

workflow GENERATE_HICUP_REPORT {
    take:
    ch_truncate_reports
    ch_map_reports
    ch_filter_reports
    ch_deduplicate_reports

    main:
    ch_truncate_reports
        .map {
             meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = meta.id - ~/(_[a-z]{4})$/
                [ meta_clone, file ]
        }
        .groupTuple ()
        .set { ch_grouped_truncate_reports }

    ch_map_reports
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = meta.id - ~/(_[a-z]{4})$/
                [ meta_clone, file ]
        }
        .groupTuple ()
        .set { ch_grouped_map_reports }

    ch_filter_reports
        .map {
            meta, summary, distribution ->
                def meta_clone = meta.clone()
                meta_clone.id = meta.id - ~/(_[a-z]{4})$/
                [ meta_clone, [ summary, distribution ] ]
        }
        .groupTuple ()
        .set { ch_grouped_filter_reports }

    ch_deduplicate_reports
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = meta.id - ~/(_[0-9a-zA-Z]*)$/
                [ meta_clone, file ]
        }
        .groupTuple ()
        .set { ch_grouped_deduplicate_reports }

    ch_grouped_truncate_reports
        .join ( ch_grouped_map_reports )
        .join ( ch_grouped_filter_reports )
        .join ( ch_grouped_deduplicate_reports )
        .map {
            meta, files ->
            def flatit = it.flatten()
            return tuple(flatit[0], flatit[1..-1])
        }
        .set{ ch_hicup_reports }

    HICUP_GENERATE_REPORT ( ch_hicup_reports )

    emit:
    multiqc = HICUP_GENERATE_REPORT.out.multiqc
}
