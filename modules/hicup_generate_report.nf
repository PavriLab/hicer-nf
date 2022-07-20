process HICUP_GENERATE_REPORT {

    tag "$meta.id"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')

    output:
    path("*html"),                   emit: html
    path("HiCUP_summary_report*"),   emit: multiqc

    script:
    def resourceDir = "${NXF_HOME}/assets/pavrilab/hicer-nf/resource"

    """
    hicupReportMerger.py \
        -o . \
        ${reports} \
        ${resourceDir}/hicup_report_template.html \
        ${meta.id}
    """
}
