process TRIM_GALORE {

    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_val_*.fq.gz"), emit: reads
    path "*fastqc.{zip,html}",                      emit: fastqc
    path "*trimming_report.txt",                    emit: reports

    shell:
    println(reads[0].name)
    '''
    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --illumina \
                --gzip \
                --output_dir !{meta.id} \
                --basename !{meta.id}_trimmed \
                --cores !{task.cpus} \
                !{reads}

    mv !{reads[0].name}_trimming_report.txt !{meta.id}_trimmed_val_1.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{meta.id}_trimmed_val_1/g' !{meta.id}_trimmed_val_1.fq.gz_trimming_report.txt
    mv !{reads[1].name}_trimming_report.txt !{meta.id}_trimmed_val_2.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{meta.id}_trimmed_val_2/g' !{meta.id}_trimmed_val_2.fq.gz_trimming_report.txt
    '''
}
