include { MAKE_BASE_MATRIX } from '../modules/cooler_make_base_matrix.nf'
include { MAKE_MCOOL       } from '../modules/cooler_make_mcool.nf'
include { BALANCE_MCOOL    } from '../modules/cooler_balance_mcool.nf'

workflow COOLER_MAKE_MATRIX {
  take:
  ch_pairs

  main:
  ch_pairs | \
  MAKE_BASE_MATRIX | \
  MAKE_MCOOL | \
  BALANCE_MCOOL

  emit:
  
}
