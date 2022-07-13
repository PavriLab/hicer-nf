include { DIGEST_GENOME       } from '../modules/hicup_digest_genome.nf'
include { BOWTIE2_BUILD_INDEX } from '../modules/bowtie2_build_index.nf'

workflow PREPARE_GENOME {
  take:
  ch_genome

  main:
  ch_genome | \
  //how do we actually fork this channel here to account for both processes

  emit:
  
}
