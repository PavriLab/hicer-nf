process TRIM_GALORE {

    tag { parameters.name }

    input:
    val(parameters) from samplesChannel

    output:
    file "${parameters.name}/*_fastqc.{zip,html}" into fastqcResults
    file "${parameters.name}/*trimming_report.txt" into trimgaloreResults
    tuple val("${parameters.name}"), file("${parameters.name}/*_trimmed_val_1.fq.gz"), file("${parameters.name}/*_trimmed_val_2.fq.gz") into resultsTrimming

    shell:
    lastPath = parameters.read1.lastIndexOf(File.separator)
    read1Base = parameters.read1.substring(lastPath+1)
    lastPath = parameters.read2.lastIndexOf(File.separator)
    read2Base = parameters.read2.substring(lastPath+1)
    trimDir = parameters.name

    """
    mkdir -p !{parameters.name}

    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --illumina \
                --gzip \
                --output_dir !{parameters.name} \
                --basename !{parameters.name}_trimmed \
                --cores !{task.cpus} \
                !{parameters.read1} \
                !{parameters.read2}

    mv !{trimDir}/!{read1Base}_trimming_report.txt !{trimDir}/!{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_1/g' !{trimDir}/!{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    mv !{trimDir}/!{read2Base}_trimming_report.txt !{trimDir}/!{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_2/g' !{trimDir}/!{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    """
}
