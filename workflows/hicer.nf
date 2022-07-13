include { PREPARE_GENOME     } from '../subworkflows/prepare_genome.nf'
include { HIC                } from '../subworkflows/hic.nf'
include { MICROC             } from '../subworkflows/microc.nf'
include { MAKE_PAIRS_FILE    } from '../modules/make_pairs_file.nf'
include { COOLER_MAKE_MATRIX } from '../subworkflows/cooler_make_matrix.nf'
include { JUICER_MAKE_MATRIX } from '../subworkflows/juicer_make_matrix.nf'
include { MULTIQC            } from '../modules/multiqc.nf'

workflow HICER {
  take:
  // check how to take more than one Channel
  ch_fastqs, ch_genome

  main:
  // again find out how to route everything

  emit:

}
