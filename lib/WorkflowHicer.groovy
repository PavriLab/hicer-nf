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
        resolutionsList = new ArrayList<Integer>()

        for (String s: resolutions_string.split(',')) {
            resolutionsList.add(s.toInteger())
        }

        resolutionsList.sort()
        resolutionsList.unique()

        sb = new StringBuilder()
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
        resolutionsList = new ArrayList<Integer>()
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
