process BOWTIE2_BUILD_INDEX {

    tag "${genome_base}"
    memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
    time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

    input:
    tuple val(genome_base), file(genomeFasta)
    val(genomeSizeType)

    output:
    tuple val(genome_base), path("bowtie2Index")

    script:
    def largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
    
    """
    mkdir bowtie2Index

    bowtie2-build \
        ${genomeFasta} \
        bowtie2Index/${genome_base} \
        --threads ${task.cpus} \
        ${largeIndexFlag}
    """
}
