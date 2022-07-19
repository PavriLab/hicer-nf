// similar to the nf-core module GUNZIP
process GUNZIP {
    tag "$archive"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$gunzip"), emit: gunzip

    script:
    gunzip = archive.toString() - '.gz'
    """
    gunzip \\
        -f \\
        $archive
    """
}
