import groovy.xml.XmlParser

class WorkflowHicer {
    // public methods
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.genome && !params.fasta) {
            log.error "Neither genome nor fasta file are specified"
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
    public static ArrayList<Integer> makeResolutionList(resolutions_string) {
        def resolutionsList = new ArrayList<Integer>()
        for (String s: resolutions_string.split(',')) {
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
              def result = parser.parse( chromSizesFile )
              result
                  .children()
                  .each{ genomeSize += it.@totalBases.toLong() }

        } else {
              new File( chromSizesFile ).eachLine{ line -> genomeSize += line.split('\t')[1].toLong() }
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
        log.info " Output Directory         : ${params.outdir}"
        log.info " ======================"
        log.info ""
    }

    public static ArrayList distributeMetaSingle(tuple) {
        return distributeMeta(tuple)
    }

    public static ArrayList distributeMetaPaired(tuple) {
        def transformed_tuple = [ tuple[0], pairFiles(tuple[1]) ]
        return distributeMeta(transformed_tuple)
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

    //
    // distributes meta to all items in a collection
    //
    private static ArrayList distributeMeta(tuple) {
        def meta    = tuple[0]
        def items   = tuple[1]
        def distributed = []
        for (item in items) {
            distributed.add( [ meta, item ] )
        }
        return distributed
    }

    //
    // paires up files in a list of files where two successive files are a pair
    //
    private static ArrayList pairFiles(files) {
        def read1 = []
        def read2 = []
        files.eachWithIndex { v, ix -> ( ix & 1 ? read2 : read1 ) << v }
        return [read1, read2].transpose()
    }

}
