process BALANCE_MCOOL {

    publishDir  path: "${params.outputDir}/${name}/matrices/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "${name}.mcool"

    tag { name }

    input:
    tuple val(name), file(mcool) from resultsZoomifyBase

    output:
    tuple val(name), file("${mcool}") into resultsMcoolNormalizer

    shell:

    '''
    balanceMultiCooler.py -m !{mcool} -p !{task.cpus}
    '''
}
