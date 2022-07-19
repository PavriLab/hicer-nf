process BOWTIE2_BUILD_INDEX {

  tag "${bwt2_base}"
  memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
  time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

  input:
  file(genomeFasta)
  val(genomeSizeType)

  output:
  file("bowtie2Index")

  script:
  largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
  lastPath = genomeFasta.lastIndexOf(File.separator)
  bwt2_base = genomeFasta.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

  """
  mkdir bowtie2Index

  bowtie2-build ${genomeFasta} \
                bowtie2Index/${bwt2_base} \
                --threads !{task.cpus} \
                ${largeIndexFlag}
  """

}
