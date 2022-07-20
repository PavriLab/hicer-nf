process HICUP_GENERATE_REPORT {

    tag "$meta.id"

    input:
    tuple val(meta), file(reports)

    output:
    path("${meta.id}/*html"),                   emit: html
    path("${meta.id}/HiCUP_summary_report*"),   emit: multiqc

    script:
    def resourceDir = "${NXF_HOME}/assets/pavrilab/hicer-nf/resource"

    """
    mkdir ${meta.id}
    hicupReportMerger.py \
        -o ${meta.id} \
        ${reports} \
        ${resourceDir}/hicup_report_template.html \
        ${meta.id}
    """
}
