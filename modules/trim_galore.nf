process TRIM_GALORE {

    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_val_*.fq.gz"), emit: reads
    path "*fastqc.{zip,html}",                      emit: fastqc
    path "*trimming_report.txt",                    emit: reports

    script:
    read1       = reads[0]
    read2       = reads[1]
    lastPath    = read1.lastIndexOf(File.separator)
    read1Base   = read1.substring(lastPath+1)
    lastPath    = read2.lastIndexOf(File.separator)
    read2Base   = read2.substring(lastPath+1)

    """
    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --illumina \
                --gzip \
                --output_dir ${meta.id} \
                --basename ${meta.id}_trimmed \
                --cores ${task.cpus} \
                ${read1} \
                ${read2}

    mv ${read1Base}_trimming_report.txt ${parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: ${meta.id}_trimmed_val_1/g' ${meta.id}_trimmed_val_1.fq.gz_trimming_report.txt
    mv ${read2Base}_trimming_report.txt ${meta.id}_trimmed_val_2.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: ${meta.id}_trimmed_val_2/g' ${meta.id}_trimmed_val_2.fq.gz_trimming_report.txt
    """
}
