process RESPLIT_FILTERED_PAIRS {

    tag { name }

    input:
    tuple val(name), file(filterSams) from resplitInputChannel

    output:
    file "${name}/${name}*.sam" into resultsResplit

    shell:
    '''
    mkdir !{name}
    resplitByChromosome.py -i !{filterSams} \
                           -o !{name}/!{name}
    '''
}
