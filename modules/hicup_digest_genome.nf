process DIGEST_GENOME {

    tag "${fasta}"

    input:
    file(genomeFasta)
    val(genomeName)
    val(re_pattern)

    output:
    path "Digest*.txt", emit: digest

    script:
    """
    hicup_digester \
        --genome ${genomeName} \
        --re1 ${re_pattern} \
        ${genomeFasta}
    """
}
