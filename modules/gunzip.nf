// similar to the nf-core module GUNZIP
process GUNZIP {
    tag "${archive}"

    input:
    file(archive)

    output:
    tuple val(genome_base), path("${gunzip}"), emit: gunzip

    script:
    def gunzip = archive.toString() - '.gz'
    def genome_base = archive.getSimpleName() - '.fa'
    """
    gunzip \
        -f \
        ${archive}
    """
}
