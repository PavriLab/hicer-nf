process MAKE_MCOOL {

    tag { name }

    input:
    tuple val(name), file(basematrix) from resultsBaseBuilder

    output:
    tuple val(name), file("${name}/${name}.mcool") into resultsZoomifyBase

    shell:
    '''
    mkdir -p !{name}

    cooler zoomify -p !{task.cpus} \
                   -r !{resolutions} \
                   -o !{name}/!{name}.mcool \
                   !{basematrix}
    '''
}
