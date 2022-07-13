include { MAKE_MATRIX } from '../modules/juicer_make_matrix.nf'

workflow JUICER_MAKE_MATRIX {
  take:
  ch_pairs

  main:
  ch_pairs | \
  MAKE_MATRIX

  emit:
  
}
