process HICUP_FILTER_PAIRS {

    tag "$meta.id"
    memory = { genomeSizeType == 'large' ? 40.GB * task.attempt : 10.GB * task.attempt }

    input:
    tuple val(meta), file(alignments)
    file(digest)

    output:
    tuple val(meta), path("${meta.id}_1_2.filt.sam"),                           emit: alignments
    tuple val(meta), path("*summary*.txt"), path("*.ditag_size_distribution"),  emit: reports

    shell:
    bin = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
    '''
    # set PERL5LIB to make hicup_module.pm available for modified hicup_filter
    LINK=$(which hicup)
    HICUPPATH=$(readlink -f $LINK)
    export PERL5LIB="$(dirname $HICUPPATH)"

    !{bin}/hicup_filter \
        --outdir . \
        --digest !{digest} \
        !{alignments}
    '''
}
