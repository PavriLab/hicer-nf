process BOWTIE2_BUILD_INDEX {

  tag "${bwt2_base}"
  memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
  time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

  input:
  file(fasta) from fastaForBowtie2

  output:
  file("bowtie2Index") into bowtie2Index

  shell:
  largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
  """
  mkdir bowtie2Index

  bowtie2-build ${fasta} \
                bowtie2Index/${bwt2_base} \
                --threads !{task.cpus} \
                ${largeIndexFlag}
  """

}
