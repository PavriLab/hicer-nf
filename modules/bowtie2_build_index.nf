process BOWTIE2_BUILD_INDEX {

    tag "${bwt2_base}"
    memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
    time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

    input:
    file(genomeFasta)
    val(genomeSizeType)

    output:
    path("bowtie2Index")

    script:
    def largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
    def bwt2_base = genomeFasta.getSimpleName()
    """
    mkdir bowtie2Index

    bowtie2-build \
        ${genomeFasta} \
        bowtie2Index/${bwt2_base} \
        --threads ${task.cpus} \
        ${largeIndexFlag}
    """
}
