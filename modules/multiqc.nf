process MULTIQC {

    tag "multiqc"

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
