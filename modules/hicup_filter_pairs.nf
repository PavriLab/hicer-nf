process HICUP_FILTER_PAIRS {

    tag { splitName }
    memory = { genomeSizeType == 'large' ? 40.GB * task.attempt : 10.GB * task.attempt }

    input:
    tuple val(splitName), file(splitSam) from resultsHicupMapper
    file(digest) from hicupDigestIndex.collect()

    output:
    file("${splitName}/${splitName}_1_2.filt.sam") into resultsHicupFilter
    tuple val(splitName), file("${splitName}/*summary*.txt"), file("${splitName}/*.ditag_size_distribution") into hicupFilterReportChannel

    shell:
    bin = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
    '''
    mkdir !{splitName}

    # set PERL5LIB to make hicup_module.pm available for modified hicup_filter
    LINK=$(which hicup)
    HICUPPATH=$(readlink -f $LINK)
    export PERL5LIB="$(dirname $HICUPPATH)"

    !{bin}/hicup_filter --outdir !{splitName} \
                        --digest !{digest} \
                        !{splitSam}
    '''
}
