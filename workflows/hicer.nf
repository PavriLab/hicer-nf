/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowHicer.initialise(params, log)

if (params.genome && params.genomes && !params.igenomes_ignore) {
    igenomes_bowtie2    = WorkflowHicer.getGenomeAttribute(params, 'bowtie2')
    igenomes_fasta      = WorkflowHicer.getGenomeAttribute(params, 'fasta')
    igenomes_chromSizes = WorkflowHicer.getGenomeAttribute(params, 'chromSizes')

} else {
    igenomes_bowtie2 = ''
    igenomes_fasta = ''
    igenomes_chromSizes = ''
}

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input,
    params.fasta,
    params.chromSizes,
    igenome_bowtie2,
    igenome_fasta,
    igenome_chromSizes
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.resolutions) {
    resolutions = WorkflowHicer.makeResolutionsUnique(params.defaultResolutions + ',' + params.resolutions)

} else {
    resolutions = params.defaultResolutions
}

resolutionsList = WorkflowHicer.makeResolutionList(resolutions)
baseResolution = resolutionsList[0]

// setting up for prepare genome subworkflow
def prepare_genome_for_tools = []

// if we do not have --genome
if (!params.genome) {
    // bowtie2 index
    if (params.fasta) {
        prepare_genome_for_tools << 'bowtie2'

    } else {
        log.error "Neither --genome nor --fasta are specified but needed for bowtie2 index."
        System.exit(1)
    }

    // HICUP genome digest
    if (params.fasta && params.chromSizes) {
        if (params.re) {
            prepare_genome_for_tools << 'hicup'

        // if digest pattern is not specified we assume Micro-C run
        } else {
            log.warn "--re not specified. Assuming you are analysing Micro-C data"
        }

    // if genome info is incomplete we exit
    } else {
        log.error "--fasta or --chromSizes not given and --genome not specified"
        System.exit(1)
    }

// if --genome is specified we check if everything is there
} else {
    if (!igenomes_bowtie2) {
        log.info "Bowtie2 index not found in igenomes config file. Computing from igenomes_fasta"
        prepare_genoe_for_tools << 'bowtie2'
    }

    prepare_genome_for_tools << 'hicup'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_GENOME     } from '../subworkflows/prepare_genome.nf'
include { HIC                } from '../subworkflows/hic.nf'
include { MICROC             } from '../subworkflows/microc.nf'
include { MAKE_PAIRS_FILE    } from '../modules/make_pairs_file.nf'
include { COOLER_MAKE_MATRIX } from '../subworkflows/cooler_make_matrix.nf'
include { JUICER_MAKE_MATRIX } from '../subworkflows/juicer_make_matrix.nf'
include { MULTIQC            } from '../modules/multiqc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow HICER {

    main:
    // again find out how to route everything

    emit:

}
