process MAKE_MATRIX {

    tag "$meta.id"

    when:
    !params.skip_juicer

    input:
    tuple val(meta), path(pairs)
    file(genomeSizes)
    val(resolutions)

    output:
    file("${meta.id}.hic")

    script:
    juicerGenomes = [
        'hg18', 'hg19', 'hg38', 'dMel',
        'mm9', 'mm10', 'anasPlat1', 'bTaurus3',
        'canFam3', 'equCab2', 'galGal4', 'Pf3D7',
        'sacCer3', 'sCerS288c', 'susScr3', 'TAIR10'
    ]
    genome = juicerGenomes.contains(params.genome) ? params.genome : genomeSizes
    juicerPath = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
    """
    mkdir -p !{name}

    java -Xmx${task.memory.toGiga()}G -jar ${juicerPath}/juicer_tools_1.22.01.jar pre \
       -r ${resolutions} \
       -k KR,GW_KR \
       --threads ${task.cpus} \
       ${pairs[0]} \
       ${meta.id}.hic \
       ${genome}
    """
}
