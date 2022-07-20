process HICUP_GENERATE_REPORT {

    tag "$meta.id"

    publishDir path: "${params.outputDir}/${meta.id}/QC/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${meta.id}/*html",
               saveAs: {
                   file ->
                        if (file.endsWith(".html")) file(file).getName()
               }

    input:
    tuple val(meta), file(reports)

    output:
    path("${meta.id}/*html"),                   emit: html
    path("${meta.id}/HiCUP_summary_report*"),   emit: multiqc

    script:
    resourceDir = "${NXF_HOME}/assets/pavrilab/hicer-nf/resource"

    """
    mkdir ${meta.id}
    hicupReportMerger.py \
        -o ${meta.id} \
        ${reports} \
        ${resourceDir}/hicup_report_template.html \
        ${meta.id}
    """
}
