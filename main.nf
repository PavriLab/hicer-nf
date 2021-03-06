#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2020 Tobias Neumann, Daniel Malzl
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

def helpMessage() {
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
        --bowtie2Index   Optional: Path to bowtie2 index
        --hicupDigest    Restriction site digest file for HICUP

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

params.help = false
if (params.genome && params.genomes && !params.igenomes_ignore) {
    igenomes_bowtie2 = params.genomes[ params.genome ].bowtie2 ?: false
    igenomes_fasta = params.genomes[ params.genome ].fasta ?: false
    igenomes_chromSizes = params.genomes[ params.genome ].chromSizes ?: false

} else {
    igenomes_bowtie2 = false
    igenomes_fasta = false
    igenomes_chromSizes = false
}

if (params.help) {
    helpMessage()
    exit 0
}

if (params.resolutions) {
  tmpResolutionsString = params.defaultResolutions + ',' + params.resolutions
  resolutionsList = new ArrayList<Integer>()

  for (String s: tmpResolutionsString.split(',')) {
    resolutionsList.add(s.toInteger())
  }

  resolutionsList.sort()
  resolutionsList.unique()
  baseResolution = resolutionsList[0]

  sb = new StringBuilder()
  for (Integer i: resolutionsList) {
    if (sb.length() > 0 ) {
      sb.append(",");
    }
    sb.append(i.toString())
  }

  resolutions = sb.toString()

} else {
  resolutions = params.defaultResolutions
  resolutionsList = new ArrayList<Integer>()

  for (String s: resolutions.split(',')) {
    resolutionsList.add(s.toInteger())
  }

  baseResolution = resolutionsList[0]

}

if (!params.bowtie2Index || !params.hicupDigest) {
  if (!params.fasta && !igenomes_fasta) {
    if (params.re) {
      exit 1, "Fasta needed for Bowtie2Index or HICUP Digest not specified!"

    } else {
      exit 1, "Fasta needed for Bowtie2Index not specified"

    }

  } else if (params.fasta) {
    if (params.re) {
      Channel
          .fromPath(params.fasta, checkIfExists: true)
          .ifEmpty { exit 1, "Fasta needed for Bowtie2Index or HICUP Digest given but not found at ${params.fasta}"}
      fastaFile = params.fasta

    } else {
      Channel
          .fromPath(params.fasta, checkIfExists: true)
          .ifEmpty { exit 1, "Fasta needed for Bowtie2Index given but not found at ${params.fasta}"}
      fastaFile = params.fasta

    }

  } else {
    if (params.re) {
      Channel
          .fromPath(igenomes_fasta, checkIfExists: true)
          .ifEmpty { exit 1, "Fasta needed for Bowtie2Index or HICUP Digest but not given and not found in igenomes or igenomes not present."}
      fastaFile = igenomes_fasta

    } else {
      Channel
          .fromPath(igenomes_fasta, checkIfExists: true)
          .ifEmpty { exit 1, "Fasta needed for Bowtie2Index but not given and not found in igenomes or igenomes not present."}
      fastaFile = igenomes_fasta

    }

  }

} else {
  if (params.re) {
    log.info "bowtie2Index and hicupDigest are specified explicitly. Fasta file will not be used if given!"
    fastaFile = "not used due to --hicupDigest and --bowtie2Index set"

  } else {
    log.info "bowtie2Index is specified explicitly. Fasta file will not be used if given!"
    fastaFile = "not used due to --bowtie2Index set"

  }

}

if (!params.bowtie2Index) {
  if (params.fasta) {
    lastPath = fastaFile.lastIndexOf(File.separator)
    bwt2_base = fastaFile.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

    fastaForBowtie2 = Channel
                          .fromPath(fastaFile)
    bowtie2IndexFile = 'computed from fasta'
    makeBowtie2Index = true

  }  else if (igenomes_bowtie2) {
    lastPath = igenomes_bowtie2.lastIndexOf(File.separator)
    bwt2_dir = igenomes_bowtie2.substring(0,lastPath+1)
    bwt2_base = igenomes_bowtie2.substring(lastPath+1)

    bowtie2Index = Channel
                      .fromPath(bwt2_dir , checkIfExists: true)
                      .ifEmpty { exit 1, "Genome index: Provided index not found: ${igenomes_bowtie2}" }
    bowtie2IndexFile = igenomes_bowtie2
    makeBowtie2Index = false

  } else {
    exit 1, "Bowtie2Index was not specified! Use either --fasta or make sure an igenomes database is properly configured."
  }

} else {
  bowtie2IndexFile = params.bowtie2Index
  lastPath = params.bowtie2Index.lastIndexOf(File.separator)
  bwt2_dir = params.bowtie2Index.substring(0,lastPath+1)
  bwt2_base = params.bowtie2Index.substring(lastPath+1)
  bowtie2Index = Channel
                    .fromPath(bwt2_dir , checkIfExists: true)
                    .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bowtie2Index}" }
  makeBowtie2Index = false

}

if (params.hicupDigest) {
  Channel
      .fromPath(params.hicupDigest, checkIfExists: true)
      .ifEmpty { exit 1, "HICUP Digest not found: ${params.hicupDigest}" }
      .set { hicupDigestIndex }
  hicupDigestFile = params.hicupDigest
  digestFasta = false

} else if (params.re) {
  fastaForHicupDigest = Channel
                            .fromPath(fastaFile)
  hicupDigestFile = 'computed from fasta'
  digestFasta = true

} else {
    log.info "HICUP digest and re motif not given. Results in assumption of micro-C protocol and skipping of truncation and alteration of filtering"
    digestFasta = false
}

if (params.chromSizes) {
  chromSizesFile = params.chromSizes

} else if (igenomes_chromSizes) {
  chromSizesFile = igenomes_chromSizes

} else {
  exit 1, "--chromSizes not specified!"
}

genomeSize = 0
if (chromSizesFile.endsWith('xml')) {
  def parser = new XmlParser()
  result = parser.parse(chromSizesFile)
  result
      .children()
      .each { genomeSize += it.@totalBases.toLong() }

  xml2tsvChannel = Channel
                      .fromPath(chromSizesFile, checkIfExists: true)
                      .ifEmpty { exit 1, "chromSize file not found at ${chromSizesFile}" }
  convertChromSizes = true

} else {
  file(chromSizesFile)
      .eachLine{ str -> genomeSize += str.split('\t')[1].toLong() }

  Channel
      .fromPath(chromSizesFile, checkIfExists: true)
      .ifEmpty { exit 1, "chromSize file not found at ${chromSizesFile}" }
      .into { chromSizeChannelPairix; chromSizeChannelCooler; chromSizeChannelJuicer }
  convertChromSizes = false

}

genomeSizeType = genomeSize > 4000000000 ? 'large' : 'small'
genomeName = params.genome ? params.genome : file(fastaFile).getSimpleName()

if (params.re) {
  log.info ""
  log.info " parameters "
  log.info " ======================"
  log.info " Samples List             : ${params.samples}"
  log.info " Resolutions              : ${resolutions}"
  log.info " baseResolution           : ${baseResolution}"
  log.info " re                       : ${params.re}"
  log.info " Genome                   : ${genomeName}"
  log.info " Genome Size              : ${genomeSizeType}"
  log.info " Fasta                    : ${fastaFile}"
  log.info " ChromSizes               : ${chromSizesFile}"
  log.info " Bowtie2 Index            : ${bowtie2IndexFile}"
  log.info " HICUP Digest             : ${hicupDigestFile}"
  log.info " Output Directory         : ${params.outputDir}"
  log.info " ======================"
  log.info ""

} else {
  log.info ""
  log.info " parameters "
  log.info " ======================"
  log.info " Samples List             : ${params.samples}"
  log.info " Resolutions              : ${resolutions}"
  log.info " baseResolution           : ${baseResolution}"
  log.info " minMapDistance           : ${params.minMapDistance}"
  log.info " Genome                   : ${genomeName}"
  log.info " Genome Size              : ${genomeSizeType}"
  log.info " Fasta                    : ${fastaFile}"
  log.info " ChromSizes               : ${chromSizesFile}"
  log.info " Bowtie2 Index            : ${bowtie2IndexFile}"
  log.info " Output Directory         : ${params.outputDir}"
  log.info " ======================"
  log.info ""

}

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .into { samplesChannel ; optionalDiscoveryChannel }

if (convertChromSizes) {
  process xml2tsv {

    tag "xml2tsv"

    input:
    file(chromSizeXML) from xml2tsvChannel

    output:
    file("chromSizes.tsv") into (chromSizeChannelPairix, chromSizeChannelCooler, chromSizeChannelJuicer)

    shell:
    '''
    xml2tsv.py !{chromSizeXML} chromSizes.tsv
    '''
  }
}


if (digestFasta) {
  process makeHicupDigest {

    tag "${fasta}"

    input:
    file(fasta) from fastaForHicupDigest

    output:
    file("Digest*.txt") into hicupDigestIndex

    shell:
    """
    echo !{task.memory}
    echo !{task.cpus}
    hicup_digester --genome !{genomeName} --re1 !{params.re} !{fasta}
    """
  }
}

if (makeBowtie2Index) {
  process buildBowtie2Index {

    tag "${bwt2_base}"
    memory = { genomeSizeType == 'large' ? 100.GB * task.attempt : 20.GB * task.attempt }
    time = { genomeSizeType == 'large' ? 8.h * task.attempt : 4.h * task.attempt }

    input:
    file(fasta) from fastaForBowtie2

    output:
    file("bowtie2Index") into bowtie2Index

    shell:
    largeIndexFlag = genomeSizeType == 'large' ? '--large-index' : ''
    """
    mkdir bowtie2Index

    bowtie2-build ${fasta} \
                  bowtie2Index/${bwt2_base} \
                  --threads !{task.cpus} \
                  ${largeIndexFlag}
    """

  }
}

process trim {

    tag { parameters.name }

    input:
    val(parameters) from samplesChannel

    output:
    file "${parameters.name}/*_fastqc.{zip,html}" into fastqcResults
    file "${parameters.name}/*trimming_report.txt" into trimgaloreResults
    tuple val("${parameters.name}"), file("${parameters.name}/*_trimmed_val_1.fq.gz"), file("${parameters.name}/*_trimmed_val_2.fq.gz") into resultsTrimming

    shell:
    lastPath = parameters.read1.lastIndexOf(File.separator)
    read1Base = parameters.read1.substring(lastPath+1)
    lastPath = parameters.read2.lastIndexOf(File.separator)
    read2Base = parameters.read2.substring(lastPath+1)
    trimDir = parameters.name

    """
    mkdir -p !{parameters.name}

    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --illumina \
                --gzip \
                --output_dir !{parameters.name} \
                --basename !{parameters.name}_trimmed \
                --cores !{task.cpus} \
                !{parameters.read1} \
                !{parameters.read2}

    mv !{trimDir}/!{read1Base}_trimming_report.txt !{trimDir}/!{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_1/g' !{trimDir}/!{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    mv !{trimDir}/!{read2Base}_trimming_report.txt !{trimDir}/!{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_2/g' !{trimDir}/!{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    """
}

process splitFastqs {

    tag { name }

    input:
    tuple val(name), file(fastq1), file(fastq2) from resultsTrimming

    output:
    file("${name}/${name}_*") into resultsSplitting

    shell:
    numberOfLinesPerSplit = params.readsPerSplit.toInteger() * 4
    '''
    mkdir !{name}

    # ampersand at the end of the line lets linux execute the commands in parallel
    zcat !{fastq1} | \
    split -l !{numberOfLinesPerSplit} -a 4 --additional-suffix _1.fq - !{name}/!{name}_ &

    zcat !{fastq2} | \
    split -l !{numberOfLinesPerSplit} -a 4 --additional-suffix _2.fq - !{name}/!{name}_
    '''
}

resultsSplitting
            .flatten()
            .map { file ->
                      def key = file.name.toString() - ~/(_[12]\.fq)?$/
                      return tuple(key, file) }
            .groupTuple()
            .set { truncaterInputChannel }


if (params.re) {
  process hicupTruncater {

      tag { splitName }

      input:
      tuple val(splitName), file(fastqSplitPairs) from truncaterInputChannel

      output:
      tuple val(splitName), file("${splitName}/${splitName}_*.trunc.fastq") into resultsHicupTruncater
      tuple val(splitName), file("${splitName}/*summary*.txt") into hicupTruncaterReportChannel

      shell:
      '''
      mkdir !{splitName}
      hicup_truncater --outdir !{splitName} \
                      --threads !{task.cpus} \
                      --re1 !{params.re} \
                      !{fastqSplitPairs[0]} \
                      !{fastqSplitPairs[1]}
      '''
  }

} else {
  process analyzeTrimmed {

      tag { splitName }

      input:
      tuple val(splitName), file(fastqSplitPairs) from truncaterInputChannel

      output:
      tuple val(splitName), file("${splitName}/${splitName}_*.trunc.fastq") into resultsHicupTruncater
      tuple val(splitName), file("${splitName}/*summary*.txt") into hicupTruncaterReportChannel

      shell:
      '''
      mkdir !{splitName}
      dummyReportGenerator.py -1 !{fastqSplitPairs[0]} \
                              -2 !{fastqSplitPairs[1]} \
                              -o !{splitName}/!{splitName}.dummy_truncater_summary.txt

      cp !{fastqSplitPairs[0]} !{splitName}/!{splitName}_1.trunc.fastq
      cp !{fastqSplitPairs[1]} !{splitName}/!{splitName}_2.trunc.fastq
      '''
  }

}

process hicupMapper {

    tag { splitName }
    memory = { genomeSizeType == 'large' ? 80.GB * task.attempt : 20.GB * task.attempt }

    input:
    tuple val(splitName), file(fastqTruncPairs) from resultsHicupTruncater
    file(index) from bowtie2Index.collect()

    output:
    tuple val(splitName), file("${splitName}/${splitName}_1_2.pair.sam") into resultsHicupMapper
    tuple val(splitName), file("${splitName}/*summary*.txt") into hicupMapperReportChannel

    shell:
    '''
    mkdir !{splitName}
    hicup_mapper --outdir !{splitName} \
                 --threads !{task.cpus} \
                 --format Sanger \
                 --index !{index}/!{bwt2_base} \
                 --bowtie2 $(which bowtie2) \
                 !{fastqTruncPairs[0]} \
                 !{fastqTruncPairs[1]}
    '''
}

if (params.re) {
  process hicupFilter {

      tag { splitName }
      memory = { genomeSizeType == 'large' ? 40.GB * task.attempt : 10.GB * task.attempt }

      input:
      tuple val(splitName), file(splitSam) from resultsHicupMapper
      file(digest) from hicupDigestIndex.collect()

      output:
      file("${splitName}/${splitName}_1_2.filt.sam") into resultsHicupFilter
      tuple val(splitName), file("${splitName}/*summary*.txt"), file("${splitName}/*.ditag_size_distribution") into hicupFilterReportChannel

      shell:
      bin = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
      '''
      mkdir !{splitName}

      # set PERL5LIB to make hicup_module.pm available for modified hicup_filter
      LINK=$(which hicup)
      HICUPPATH=$(readlink -f $LINK)
      export PERL5LIB="$(dirname $HICUPPATH)"

      !{bin}/hicup_filter --outdir !{splitName} \
                          --digest !{digest} \
                          !{splitSam}
      '''
  }

} else {
  process sizeFilter {

      tag { splitName }

      input:
      tuple val(splitName), file(splitSam) from resultsHicupMapper

      output:
      file("${splitName}/${splitName}_1_2.filt.sam") into resultsHicupFilter
      tuple val(splitName), file("${splitName}/*summary*.txt"), file("${splitName}/*.ditag_size_distribution") into hicupFilterReportChannel

      shell:
      '''
      mkdir !{splitName}
      filterBySize.py -i !{splitSam} \
                      --minDistance !{params.minMapDistance} \
                      -o !{splitName}/!{splitName}_1_2.filt.sam
      '''
  }

}

resultsHicupFilter
              .map { file ->
                        def key = file.name.toString() - ~/(_[a-z]{4}_1_2\.filt\.sam)?$/
                        return tuple(key, file) }
              .groupTuple()
              .set { resplitInputChannel }

process resplitFiltered {

    tag { name }

    input:
    tuple val(name), file(filterSams) from resplitInputChannel

    output:
    file "${name}/${name}*.sam" into resultsResplit

    shell:
    '''
    mkdir !{name}
    resplitByChromosome.py -i !{filterSams} \
                           -o !{name}/!{name}
    '''
}

process hicupDeduplicator {

    tag { resplitName }

    input:
    file sam from resultsResplit.flatten()

    output:
    file("${resplitName}/${resplitName}_1_2.dedup.sam") into resultsHicupDeduplicate
    tuple val(resplitName), file("${resplitName}/*summary*.txt") into hicupDeduplicatorReportChannel

    shell:
    resplitName = (sam.getName() - ~/(_1_2\.filt\.sam)?$/)

    '''
    mkdir !{resplitName}
    hicup_deduplicator --outdir !{resplitName} \
                       !{sam}

    if [ ! !{params.re} ]
    then
        getCisFractions.py -i !{resplitName}/!{resplitName}_1_2.dedup.sam \
                           -r !{resplitName}/*summary*.txt
    fi
    '''
}

resultsHicupDeduplicate
              .map { file ->
                        def key = file.name.toString() - ~/(_[0-9a-zA-Z]*_1_2\.dedup\.sam)?$/
                        return tuple(key, file) }
              .groupTuple()
              .set { catSamInputChannel }

process catSam {

    tag { name }

    input:
    tuple val(name), file(dedupSams) from catSamInputChannel

    output:
    tuple val(name), file("${name}_1_2.dedup.sam") into resultsHicup, sam2bamChannel

    shell:
    '''
    samtools view -H !{dedupSams[0]} > !{name}_1_2.dedup.sam
    cat !{dedupSams} | grep -v '^@' >> !{name}_1_2.dedup.sam
    '''
}

hicupTruncaterReportChannel
                .map{ it ->
                        def key = it[0] - ~/(_[a-z]{4})$/
                        return tuple(key, it[1]) }
                .groupTuple()
                .set { hicupTruncaterGroupedChannel }

hicupMapperReportChannel
                .map{ it ->
                        def key = it[0] - ~/(_[a-z]{4})$/
                        return tuple(key, it[1]) }
                .groupTuple()
                .set { hicupMapperGroupedChannel }

hicupFilterReportChannel
                .map{ it ->
                        def key = it[0] - ~/(_[a-z]{4})$/
                        return tuple(key, it[1], it[2]) }
                .groupTuple()
                .set { hicupFilterGroupedChannel }

hicupDeduplicatorReportChannel
                .map{ it ->
                        def key = it[0] - ~/(_[0-9a-zA-Z]*)$/
                        return tuple(key, it[1]) }
                .groupTuple()
                .set { hicupDeduplicatorGroupedChannel }

hicupTruncaterGroupedChannel
                .join(hicupMapperGroupedChannel)
                .join(hicupFilterGroupedChannel)
                .join(hicupDeduplicatorGroupedChannel)
                .map { it ->
                          def flatit = it.flatten()
                          return tuple(flatit[0], flatit[1..-1])}
                .set{ hicupReporterInputChannel }

process hicupReporter {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/QC/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*html",
               saveAs: { filename ->
                             if (filename.endsWith(".html")) file(filename).getName()
                       }

    input:
    tuple val(name), file(hicupReportFiles) from hicupReporterInputChannel

    output:
    file("${name}/*html") into htmlHicup
    file("${name}/HiCUP_summary_report*") into multiqcHicup

    shell:
    resourceDir = "${NXF_HOME}/assets/pavrilab/hicer-nf/resource"
    '''
    mkdir !{name}
    hicupReportMerger.py -o !{name} \
                         !{hicupReportFiles} \
                         !{resourceDir}/hicup_report_template.html \
                         !{name}
    '''
}

process sam2bamConverter {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/bam/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*.bam",
               saveAs: { filename ->
                             if (filename.endsWith(".bam")) file(filename).getName()
                       }

    input:
    tuple val(name), file(sam) from sam2bamChannel

    output:
    file("${name}/${name}.hicup.bam") into sam2bamResults

    shell:
    '''
    mkdir -p !{name}

    samtools view -bh -@ !{task.cpus} !{sam} > !{name}/!{name}.hicup.bam
    '''
}

process pairixMaker {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/pairs/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "${name}/*pairs.gz*",
               saveAs: { filename ->
                             if (filename.endsWith(".pairs.gz")) file(filename).getName()
                             else if (filename.endsWith(".pairs.gz.px2")) file(filename).getName()
                       }

    input:
    tuple val(name), file(sam) from resultsHicup
    file(chromSizeFile) from chromSizeChannelPairix.collect()

    output:
    tuple val(name), file("${name}/${name}.pairs.gz"), file("${name}/${name}.pairs.gz.px2") into resultsPairixBaseBuilder, resultsPairixJuicer

    shell:
    '''
    mkdir -p !{name}

    sam2pairs.py !{sam} | \
        paste - - | \
        awk 'BEGIN{ FS = "\t"; OFS = "\t" }{ print $1,$2,$3,$6,$7,$4,$8 }' > \
        !{name}/!{name}.pairs.tmp

    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V !{chromSizeFile} > chromSizes.sort.tsv

    cooler csort -c1 2 -c2 4 \
                 -p1 3 -p2 5 \
                 -p !{task.cpus} \
                 !{name}/!{name}.pairs.tmp \
                 chromSizes.sort.tsv

    # add generic header to make pairix compatible with juicer prefix
    echo "## pairs format v1.0" > !{name}/!{name}.pairs
    echo "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2" >> !{name}/!{name}.pairs
    zcat !{name}/!{name}.pairs.tmp.blksrt.gz >> !{name}/!{name}.pairs

    bgzip !{name}/!{name}.pairs
    pairix -p pairs !{name}/!{name}.pairs.gz
    '''
}

process juicerHic {

  tag { name }

  publishDir path: "${params.outputDir}/${name}/matrices/",
             mode: 'copy',
             overwrite: 'true',
             pattern: "${name}/*.hic",
             saveAs: { filename ->
                           if (filename.endsWith(".hic")) file(filename).getName()
                     }
  when:
  !params.skip_juicer
  
  input:
  tuple val(name), file(pairs), file(pairsIndex) from resultsPairixJuicer
  file(chromSizeFile) from chromSizeChannelJuicer.collect()

  output:
  file("${name}/${name}.hic") into resultsJuicer

  shell:
  juicerGenomes = ['hg18', 'hg19', 'hg38', 'dMel',
                   'mm9', 'mm10', 'anasPlat1', 'bTaurus3',
                   'canFam3', 'equCab2', 'galGal4', 'Pf3D7',
                   'sacCer3', 'sCerS288c', 'susScr3', 'TAIR10']
  genome = juicerGenomes.contains(params.genome) ? params.genome : chromSizeFile
  juicerPath = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
  '''
  mkdir -p !{name}

  java -Xmx!{task.memory.toGiga()}G -jar !{juicerPath}/juicer_tools_1.22.01.jar pre \
       -r !{resolutions} \
       -k KR,GW_KR \
       --threads !{task.cpus} \
       !{pairs} \
       !{name}/!{name}.hic \
       !{genome}
  '''
}

process baseBuilder {

    tag { name }

    input:
    tuple val(name), file(pairs), file(pairsIndex) from resultsPairixBaseBuilder
    file(chromSizeFile) from chromSizeChannelCooler.collect()

    output:
    tuple val(name), file("${name}/${name}_base.cool") into resultsBaseBuilder

    shell:
    '''
    mkdir -p !{name}

    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V !{chromSizeFile} > chromSizes.sort.tsv

    cooler cload pairs --assembly !{genomeName} \
                       -c1 2 -p1 3 -c2 4 -p2 5 \
                       chromSizes.sort.tsv:!{baseResolution} \
                       !{pairs} \
                       !{name}/!{name}_base.cool
    '''

}

process zoomifyBase {

    tag { name }

    input:
    tuple val(name), file(basematrix) from resultsBaseBuilder

    output:
    tuple val(name), file("${name}/${name}.mcool") into resultsZoomifyBase

    shell:
    '''
    mkdir -p !{name}

    cooler zoomify -p !{task.cpus} \
                   -r !{resolutions} \
                   -o !{name}/!{name}.mcool \
                   !{basematrix}
    '''
}

process mcoolNormalizer {

    publishDir  path: "${params.outputDir}/${name}/matrices/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "${name}.mcool"

    tag { name }

    input:
    tuple val(name), file(mcool) from resultsZoomifyBase

    output:
    tuple val(name), file("${mcool}") into resultsMcoolNormalizer

    shell:

    '''
    balanceMultiCooler.py -m !{mcool} -p !{task.cpus}
    '''
}

process multiqc {

    tag { 'all' }

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file (fastqc: 'fastqc/*') from fastqcResults.collect()
    file (trim: 'trim/*') from trimgaloreResults.collect()
    file (hicpt: 'hicup/*') from multiqcHicup.collect()

    output:
    file "*multiqc_report.html" into multiqc_report

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -f -x *.run .
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
