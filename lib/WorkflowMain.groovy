class WorkflowMain {
    //
    // Print help to screen if required
    //
    public static void help(log) {
        log.info"""
        ================================================================
         hicer-nf
        ================================================================
         DESCRIPTION

         Basic processing of HiC data.

         Usage:
         nextflow run pavrilab/hicer-nf

         Options:
            --samples        Tab-delimited text file specifying the samples
                             to be processed. (default: 'samples.txt')
                             The following columns are required:
                                - name: name of sample
                                - read1: Read file with first read mates (R1) in fastq(.gz) format
                                - read2: Read file with second read mates (R2) in fastq(.gz) format

            --resolutions    comma-separated list of resolutions in bp to compute in addition to the default resolutions
    			                   default resolutions are 5000,10000,25000,50000,100000,250000,500000,1000000 and resolutions
    			                   specified via this parameter will be added to this list

            --re             regular expression to use for in-silico digestion by HICUP (e.g. ^GATC,MboI)
                             if not given the pipeline assumes the micro-C protocol was used and skips the truncation step
                             and filters the reads only by their minimum mapping distance

            --outputDir      Directory name to save results to. (default: 'results')

            --readsPerSplit  specifies the number of read per fastq split for parallel processing (default: 25.000.000)

            --minMapDistance minimum mapping distance between to reads of a pair (default: 500)
                             is only used in case of micro-C (i.e. --re is not given)

            --skip_juicer    if set, skips hic file generation with juicer

         References:
            --genome         Name of reference (hg38, mm10, ...)
            --fasta          Alternatively, path to genome fasta file which will be digested
            --chromSizes     tab-separated file containing chromosome names and their sizes

         Profiles:
            standard         local execution
            singularity      local execution with singularity
            cbe              CBE cluster execution with singularity

         Docker:
         zuberlab/hicer-nf:latest

         Authors:
         Tobias Neumann (tobias.neumann@imp.ac.at)
         Daniel Malzl (daniel.malzl@imp.ac.at)
        """.stripIndent()
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(log)
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        checkConfigProvided(workflow, log)

        // Check input has been provided
        if (!params.input) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }
    }

    //
    //  Warn if a -profile or Nextflow config has not been provided to run the pipeline
    //
    public static void checkConfigProvided(workflow, log) {
        if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
            log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
                "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
                "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
                "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
                "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
                "Please refer to the quick start section and usage docs for the pipeline.\n "
        }
    }
}
