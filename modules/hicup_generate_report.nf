process HICUP_GENERATE_REPORT {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/QC/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*html",
               saveAs: { filename ->
                             if (filename.endsWith(".html")) file(filename).getName()
                       }

    input:
    tuple val(name), file(hicupReportFiles) from hicupReporterInputChannel

    output:
    file("${name}/*html") into htmlHicup
    file("${name}/HiCUP_summary_report*") into multiqcHicup

    shell:
    resourceDir = "${NXF_HOME}/assets/pavrilab/hicer-nf/resource"
    '''
    mkdir !{name}
    hicupReportMerger.py -o !{name} \
                         !{hicupReportFiles} \
                         !{resourceDir}/hicup_report_template.html \
                         !{name}
    '''
}
