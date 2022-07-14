include { DIGEST_GENOME       } from '../modules/hicup_digest_genome.nf'
include { BOWTIE2_BUILD_INDEX } from '../modules/bowtie2_build_index.nf'

workflow PREPARE_GENOME {
  take:
  prepare_genome_for_tools
  genome_fasta
  re_pattern

  main:

  if ('bowtie2' in perpare_genome_for_tools) {
      BOWTIE2_BUILD_INDEX( genome_fasta )
  }

  if ('hicup' in prepare_genome_for_tools) {
      DIGEST_GENOME( genome_fasta, re_pattern)
  }

  emit:

}
