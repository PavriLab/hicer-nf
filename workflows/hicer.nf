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
    igenomes_fasta,
    igenomes_chromSizes
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.resolutions) {
    resolutions = WorkflowHicer.makeResolutionsUnique(
        params.defaultResolutions + ',' + params.resolutions
    )

} else {
    resolutions = params.defaultResolutions
}

resolutionsList = WorkflowHicer.makeResolutionList( resolutions )
baseResolution = resolutionsList[0]

// setting up for prepare genome subworkflow
def prepare_genome_for_tools = []

// if we do not have --genome
if (!params.genome) {
    // bowtie2 index
    if (params.fasta) {
        prepare_genome_for_tools << "bowtie2"

    } else {
        log.error "Neither --genome nor --fasta are specified but needed for bowtie2 index."
        System.exit(1)
    }

    // HICUP genome digest
    if (params.fasta && params.chromSizes) {
        if (params.re) {
            prepare_genome_for_tools << "hicup"

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
        prepare_genome_for_tools << "bowtie2"
    }

    if (igenomes_chromSizes.endsWith("xml")) {
        prepare_genome_for_tools << "chromSizes"
    }

    if (params.re) {
        prepare_genome_for_tools << "hicup"
    }
}

dynamic_params = [:]
dynamic_params.genomeFasta      = params.genome ? igenomes_fasta : params.fasta
dynamic_params.genomeSizes      = params.genome ? igenomes_chromSizes : params.chromSizes
dynamic_params.bowtie2Index     = igenomes_bowtie2 ? igenomes_bowtie2 : "computed from fasta"
dynamic_params.genomeSizeType   = WorkflowHicer.getGenomeSizeType( dynamic_params.genomeSizes )
dynamic_params.genomeName       = params.genome ? params.genome : file(fastaFile).getSimpleName()
dynamic_params.baseResolution   = baseResolution
dynamic_params.resolutions      = resolutions
dynamic_params.re               = params.re

WorkflowHicer.paramsSummaryLog( params, dynamic_params, log )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK        } from '../subworkflows/input_check.nf'
include { PREPARE_GENOME     } from '../subworkflows/prepare_genome.nf'
include { CAT_FASTQ          } from '../modules/cat_fastq.nf'
include { TRIM_GALORE        } from '../modules/trim_galore.nf'
include { SPLIT_FASTQ        } from '../modules/split_fastq.nf'
include { HIC                } from '../subworkflows/hic.nf'
include { MICROC             } from '../subworkflows/microc.nf'
include { SAM_TO_BAM         } from '../modules/sam_to_bam.nf'
include { MAKE_PAIRS_FILE    } from '../modules/make_pairs_file.nf'
include { COOLER_MAKE_MATRIX } from '../subworkflows/cooler_make_matrix.nf'
include { JUICER_MAKE_MATRIX } from '../subworkflows/juicer_make_matrix.nf'
include { MULTIQC            } from '../modules/multiqc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HICER {
    ch_input = file( params.input )

    // check input sample sheet
    // adapted from nf-core/rnaseq
    INPUT_CHECK ( ch_input )
        .reads
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    // prepare genome files
    if (!prepare_genome_for_tools.isEmpty()) {
        ch_genome = PREPARE_GENOME (
            prepare_genome_for_tools,
            dynamic_params
        )

    } else {
        ch_genome = [:]
        ch_genome.index     = file( dynamic_params.bowtie2Index ).getParent()
        ch_genome.digest    = ''
        ch_genome.sizes     = file( dynamic_params.genomeSizes )
    }

    // concatenate fastqs of samples with multiple readfiles
    CAT_FASTQ ( ch_fastq.multiple )
        .reads
        .mix ( ch_fastq.single )
        .set { ch_cat_fastq }

    // read QC
    TRIM_GALORE ( ch_cat_fastq )
        .reads
        .set { ch_trim_fastq }

    // splitting fastqs for HICUP parallelization
    // flatten().collate() relies on the fact that each file is interspersed
    // by a meta object after the first map invocation
    SPLIT_FASTQ ( ch_trim_fastq )
        .reads
        .map { WorkflowHicer.distributeMetaSingle( it ) }
        .flatten ()
        .collate ( 2 )
        .map {
            meta, file ->
                def meta_clone = meta.clone()
                meta_clone.id = file.name.toString() - ~/(_[12]\.fq)?$/
                [ meta_clone, file ]
        }
        .groupTuple (by: [0])
        .set { ch_split_fastq }

    if (params.re) {
        ch_hicup = HIC (
            ch_split_fastq,
            ch_genome.index,
            ch_genome.digest,
            dynamic_params.genomeSizeType
        )

    } else {
        ch_hicup = MICROC (
            ch_split_fastq,
            ch_genome.index,
            dynamic_params.genomeSizeType
        )
    }

    SAM_TO_BAM ( ch_hicup.alignments )

    MAKE_PAIRS_FILE (
        ch_hicup.alignments,
        ch_genome.sizes
    )

    COOLER_MAKE_MATRIX (
        MAKE_PAIRS_FILE.out.pairs,
        dynamic_params.genomeName,
        dynamic_params.baseResolution,
        dynamic_params.resolutions,
        ch_genome.sizes
    )

    if (!params.skip_juicer) {
        JUICER_MAKE_MATRIX (
            MAKE_PAIRS_FILE.out.pairs,
            ch_genome.sizes,
            dynamic_params.resolutions
        )
    }

    MULTIQC (
        TRIM_GALORE.out.fastqc.collect(),
        TRIM_GALORE.out.reports.collect(),
        ch_hicup.multiqc.collect()
    )
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
