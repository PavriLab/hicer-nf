process SAM_TO_BAM {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/bam/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*.bam",
               saveAs: { filename ->
                             if (filename.endsWith(".bam")) file(filename).getName()
                       }

    input:
    tuple val(name), file(sam) from sam2bamChannel

    output:
    file("${name}/${name}.hicup.bam") into sam2bamResults

    shell:
    '''
    mkdir -p !{name}

    samtools view -bh -@ !{task.cpus} !{sam} > !{name}/!{name}.hicup.bam
    '''
}
