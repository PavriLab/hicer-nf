// similar to the nf-core module GUNZIP
process GUNZIP {
    tag "$archive"

    input:
    path(archive)

    output:
    path("$gunzip"), emit: gunzip

    script:
    gunzip = archive.toString() - '.gz'
    """
    gunzip \\
        -f \\
        $archive
    """
}
