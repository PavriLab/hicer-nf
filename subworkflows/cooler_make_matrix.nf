include { MAKE_BASE_MATRIX } from '../modules/cooler_make_base_matrix.nf'
include { MAKE_MCOOL       } from '../modules/cooler_make_mcool.nf'
include { BALANCE_MCOOL    } from '../modules/cooler_balance_mcool.nf'

workflow COOLER_MAKE_MATRIX {
    take:
    ch_pairs
    genomeName
    baseResolution
    resolutions
    genomeSizes

    main:
    MAKE_BASE_MATRIX (
        ch_pairs,
        genomeSizes,
        genomeName,
        baseResolution
    )

    MAKE_MCOOL (
        MAKE_BASE_MATRIX.out.matrix,
        resolutions
    )

    BALANCE_MCOOL (
        MAKE_MCOOL.out.matrix
    )
}
