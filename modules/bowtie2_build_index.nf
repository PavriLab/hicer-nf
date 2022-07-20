process BOWTIE2_BUILD_INDEX {

    tag "${bwt2_base}"
    memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
    time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

    input:
    file(genomeFasta)
    val(genomeSizeType)

    output:
    tuple val(bwt2_base), path("bowtie2Index"), emit: index

    script:
    largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
    bwt2_base = genomeFasta.getSimpleName()

    """
    mkdir bowtie2Index

    bowtie2-build \
        ${genomeFasta} \
        bowtie2Index/${bwt2_base} \
        --threads ${task.cpus} \
        ${largeIndexFlag}
    """
}
