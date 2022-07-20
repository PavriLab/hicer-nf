include { MAKE_MATRIX } from '../modules/juicer_make_matrix.nf'

workflow JUICER_MAKE_MATRIX {
    take:
    ch_pairs
    chromSizeFile
    resolutions

    main:
    MAKE_MATRIX (
        ch_pairs,
        chromSizeFile,
        resolutions
    )
}
