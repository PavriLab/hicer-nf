class WorkflowHicer {
    // public methods
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.bowtie2Index && !params.genome && !params.fasta) {
            log.error "Neither genome nor fasta file nor bowtie2 index are specified but needed"
            System.exit(1)
        }

        if (!params.hicupDigest && !params.genome && !params.fasta && params.re) {
            log.error "Neither genome nor fasta file nor hicup digest are specified but RE pattern"
            System.exit(1)
        }
    }
    //
    // include specified resolutions in resolution list sort them and remove duplicates
    //
    public static String makeResolutionsUnique(resolutions_string) {
        def resolutionsList = new ArrayList<Integer>()

        for (String s: resolutions_string.split(',')) {
            resolutionsList.add(s.toInteger())
        }

        resolutionsList.sort()
        resolutionsList.unique()

        def sb = new StringBuilder()
        for (Integer i: resolutionsList) {
            if (sb.length() > 0) {
                sb.append(",");
            }
            sb.append(i.toString())
        }

        return sb.toString()
    }

    //
    // convert resolutions string to list
    //
    public static ArrayList<Integer> makeResolutionList(resolution_string) {
        def resolutionsList = new ArrayList<Integer>()
        for (String s: resolutions.split(',')) {
            resolutionsList.add(s.toInteger())
        }

        return resolutionsList
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }

    //
    // Compute genome length to determine resource needs
    // larger genomes usually need more resources
    //
    public static String getGenomeSizeType(chromSizesFile) {
        def genomeSize = 0
        if (chromSizesFile.endsWith('xml')) {
              def parser = new XmlParser()
              result = parser.parse( chromSizesFile )
              result
                  .children()
                  .each{ genomeSize += it.@totalBases.toLong() }

        } else {
              file( chromSizesFile ).eachLine{ str -> genomeSize += str.split('\t')[1].toLong() }
        }

        def genomeSizeType = genomeSize > 4000000000 ? "large" : "small"
        return genomeSizeType
    }

    //
    // Print parameter summary log to screen
    //
    public static void paramsSummaryLog(params, dynamic, log) {
        log.info ""
        log.info " parameters "
        log.info " ======================"
        log.info " Samples List             : ${params.input}"
        log.info " Resolutions              : ${dynamic.resolutions}"
        log.info " baseResolution           : ${dynamic.baseResolution}"
        if (params.re) {
            log.info " re                       : ${params.re}"

        } else {
            log.info " minMapDistance           : ${params.minMapDistance}"
        }

        log.info " Genome                   : ${dynamic.genomeName}"
        log.info " Genome Size              : ${dynamic.genomeSizeType}"
        log.info " Fasta                    : ${dynamic.genomeFasta}"
        log.info " ChromSizes               : ${dynamic.genomeSizes}"
        log.info " Bowtie2 Index            : ${dynamic.bowtie2Index}"
        log.info " Output Directory         : ${params.outputDir}"
        log.info " ======================"
        log.info ""
    }

    public static ArrayList distributeMeta(tuple) {
        def meta    = tuple[0]
        def files   = tuple[1]
        def distributed = []
        for (file in files) {
            distributed.add( [ meta, file ] )
        }
        return distributed
    }

    // private methods
    //
    // check if specified genome exists
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}
