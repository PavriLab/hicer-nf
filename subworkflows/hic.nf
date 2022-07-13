include { TRIM_GALORE             } from '../modules/trim_galore.nf'
include { SPLIT_FASTQS            } from '../modules/split_fastqs.nf'
include { HICUP_TRUNCATE_READS    } from '../modules/hicup_truncate_reads.nf'
include { HICUP_MAP_READS         } from '../modules/hicup_map_reads.nf'
include { HICUP_FILTER_PAIRS      } from '../modules/hicup_filter_pairs.nf'
include { RESPLIT_FILTERED_PAIRS  } from '../modules/resplit_filtered_pairs.nf'
include { HICUP_DEDUPLICATE_PAIRS } from '../modules/hicup_deduplicate_pairs.nf'
include { HICUP_GENERATE_REPORT   } from '../modules/hicup_generate_report.nf'

workflow HIC {
  take:
  ch_fastqs

  main:
  ch_fastqs | \
  TRIM_GALORE | \
  SPLIT_FASTQS | \
  HICUP_TRUNCATE_READS | \
  HICUP_MAP_READS | \
  HICUP_FILTER_PAIRS | \
  RESPLIT_FILTERED_PAIRS | \
  HICUP_DEDUPLICATE_PAIRS | \
  HICUP_GENERATE_REPORT

  emit:

}
