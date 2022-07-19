process MULTIQC {

    tag "multiqc"

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file (fastqc: 'fastqc/*')
    file (trim: 'trim/*')
    file (hicpt: 'hicup/*')

    output:
    file "*multiqc_report.html"

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -f -x *.run .
    """
}
