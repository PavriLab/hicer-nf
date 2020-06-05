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
     nextflow run t-neumann/hicer-nf

     Options:
        --samples        Tab-delimited text file specifying the samples
                         to be processed. (default: 'samples.txt')
                         The following columns are required:
                            - name: name of sample
                            - read1: Read file with first read mates (R1) in fastq(.gz) format
                            - read2: Read file with second read mates (R2) in fastq(.gz) format

        --resolutions    comma-separated list of resolutions in bp to compute in addition to the default resolutions

        --re1            regular expression to use for in-silico digestion by HICUP (e.g. ^GATC,MboI)
        --re2            second regular expression to use for in-silico digestion in case of double digestion protocol

        --outputDir      Directory name to save results to. (Defaults to
                         'results')

     References:
        --genome         Name of reference (hg38, mm10, ...)
        --fasta          Alternatively, path to genome fasta file which will be digested
        --chromSizes     tab-separated file containing chromosome names and their sizes
        --bowtie2        Optional: Path to bowtie2 index
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
igenomes_bowtie2 = params.genomes[ params.genome ].bowtie2 ?: false
igenomes_fasta = params.genomes[ params.genome ].fasta ?: false
igenomes_chromSizes = params.genomes[ params.genome ].chromSizes ?: false

if (params.help) {
    helpMessage()
    exit 0
}

if (params.resolutions) {
    resolutions = params.defaultResolutions + "," + params.resolutions
} else {
    resolutions = params.resolutions
}

if (!params.bowtie2 || !params.hicupDigest) {
  if (!params.fasta && !igenomes_fasta) {
    exit 1, "Fasta needed for Bowtie2Index or HICUP Digest not specified!"

  } else if (params.fasta) {
    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta needed for Bowtie2Index or HICUP Digest but not found at ${params.fasta}"}
    fastaFile = params.fasta

  } else {
    Channel
        .fromPath(igenomes_fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta needed for Bowtie2Index or HICUP Digest but not given and not found at ${params.fasta}"}
    fastaFile = igenomes_fasta

  }
}

if (!params.bowtie2) {
  if (igenomes_bowtie2) {
    lastPath = igenomes_bowtie2.lastIndexOf(File.separator)
    bwt2_dir = igenomes_bowtie2.substring(0,lastPath+1)
    bwt2_base = igenomes_bowtie2.substring(lastPath+1)

    bowtie2Index = Channel
                      .fromPath(bwt2_dir , checkIfExists: true)
                      .ifEmpty { exit 1, "Genome index: Provided index not found: ${igenomes_bowtie2}" }
    bowtie2IndexFile = igenomes_bowtie2
    makeBowtie2Index = false

  }  else {
    lastPath = fastaFile.lastIndexOf(File.separator)
    bwt2_base = fastaFile.substring(lastPath+1)

    fastaForBowtie2 = Channel
                          .fromPath(fastaFile)
    bowtie2IndexFile = 'computed from fasta'
    makeBowtie2Index = true
  }
} else {
  bowtie2IndexFile = igenomes_bowtie2
  makeBowtie2Index = false

}

if (params.hicupDigest) {
  Channel
      .fromPath(params.hicupDigest, checkIfExists: true)
      .ifEmpty { exit 1, "HICUP Digest not found: ${params.hicupDigest}" }
      .into{hicupDigestIndex}
  hicupDigestFile = params.hicupDigest
  digestFasts = false

} else if (params.re1) {
  fastaForHicupDigest = Channel
                            .fromPath(fastaFile)
  hicupDigestFile = 'computed from fasta'
  digestFasta = true

} else {
    exit 1, "HICUP digest file does not exist and --re1 is not set!"

}

if (params.chromSizes) {
  chromSizesFile = params.chromSizes

} else if (igenomes_chromSizes) {
  chromSizesFile = igenomes_chromSizes

} else {
  exit 1, "--chromSizes not specified!"
}

if (chromSizesFile.endsWith('xml')) {
  xml2tsvChannel = Channel
                      .fromPath(chromSizesFile, checkIfExists: true)
                      .ifEmpty { exit 1, "chromSize file not found at ${chromSizesFile}" }
  convertChromSizes = true

} else {
  Channel
      .fromPath(chromSizesFile, checkIfExists: true)
      .ifEmpty { exit 1, "chromSize file not found at ${chromSizesFile}" }
      .into { chromSizeChannelPairix; chromSizeChannelCooler; chromSizeChannelJuicer }
  convertChromSizes = false

}

log.info ""
log.info " parameters "
log.info " ======================"
log.info " Samples List             : ${params.samples}"
log.info " Resolutions              : ${resolutions}"
log.info " re1                      : ${params.re1}"
log.info " re2                      : ${params.re2}"
log.info " Genome                   : ${params.genome}"
log.info " Fasta                    : ${fastaFile}"
log.info " ChromSizes               : ${chromSizesFile}"
log.info " Bowtie2 Index            : ${bowtie2IndexFile}"
log.info " HICUP Digest             : ${hicupDigestFile}"
log.info " Output Directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

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
    re2 = params.re2 ? "--re2 ${params.re2}": ''
    """
    hicup_digester --genome !{params.genome} --re1 !{params.re1} !{re2} !{fasta}
    """
  }
}

if (makeBowtie2Index) {
  process buildBowtie2Index {

    tag "${bwt2_base}"

    input:
    file(fasta) from fastaForBowtie2

    output:
    file("bowtie2Index") into bowtie2Index

    shell:
    bwt2_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir bowtie2Index

    bowtie2-build ${fasta} bowtie2Index/${bwt2_base} --threads !{task.cpus}
    """

  }
}

process trim {

    tag { parameters.name }

    input:
    val(parameters) from samplesChannel

    output:
    file "*_fastqc.{zip,html}" into fastqcResults
    file "*trimming_report.txt" into trimgaloreResults
    tuple val("${parameters.name}"), file('*_trimmed_val_1.fq.gz'), file('*_trimmed_val_2.fq.gz') into resultsTrimming

    shell:
    lastPath = parameters.read1.lastIndexOf(File.separator)
    read1Base = parameters.read1.substring(lastPath+1)
    lastPath = parameters.read2.lastIndexOf(File.separator)
    read2Base = parameters.read2.substring(lastPath+1)

    """
    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --illumina \
                --gzip \
                --basename !{parameters.name}_trimmed \
                --cores !{task.cpus} \
                !{parameters.read1} \
                !{parameters.read2}

    mv !{read1Base}_trimming_report.txt !{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_1/g' !{parameters.name}_trimmed_val_1.fq.gz_trimming_report.txt
    mv !{read2Base}_trimming_report.txt !{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{parameters.name}_trimmed_val_2/g' !{parameters.name}_trimmed_val_2.fq.gz_trimming_report.txt
    """
}

process hicup {

    tag { name }

    publishDir  path: "${params.outputDir}/QC/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "*/*html"

    input:
    file index from bowtie2Index.collect()
    file digest from hicupDigestIndex.collect()
    tuple val(name), file(fastq1), file(fastq2) from resultsTrimming

    output:
    tuple val(name), file("${name}/*sam") into resultsHicup
    file("${name}/*html") into htmlHicup
    file("${name}/HiCUP_summary_report*") into multiqcHicup

    shell:

    '''
    mkdir -p !{name}

    hicup --bowtie2 $(which bowtie2) \
          --index !{index}/!{bwt2_base} \
          --digest !{digest} \
          --format Sanger \
          --outdir !{name} \
          --threads !{task.cpus} \
          !{fastq1} \
          !{fastq2}

    mv !{name}/*sam !{name}/!{name}.hicup.sam

    sed -i 's/^.*sam\t/!{name}.hicup.sam\t/g' !{name}/HiCUP_summary_report*txt

    mv !{name}/HiCUP_summary_report*txt !{name}/HiCUP_summary_report_!{name}.txt

    sed -i 's/HiCUP Processing Report - [^<]*/HiCUP Processing Report - !{name}/g' !{name}/*.HiCUP_summary_report.html
    sed -i 's/WRAP CHAR>[^<]*/WRAP CHAR>!{name}/g' !{name}/*.HiCUP_summary_report.html

    '''
}

process pairixMaker {

    tag { name }

    publishDir  path: "${params.outputDir}/${name}/pairs/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "*/*pairs.gz*"

    input:
    tuple val(name), file(sam) from resultsHicup
    file(chromSizeFile) from chromSizeChannelPairix

    output:
    tuple val(name), file("${name}/${name}.pairs.gz"), file("${name}/${name}.pairs.gz.px2") into resultsPairixCooler, resultsPairixJuicer


    shell:
    '''
    mkdir -p !{name}

    samtools view !{sam} | \
        awk 'BEGIN{ FS = "\t"; OFS = "\t" }{ print $1,$3,$4,and($2, 16)?"-":"+"; }' | \
        paste - - | \
        awk 'BEGIN{ FS = "\t"; OFS = "\t" }{ print $1,$2,$3,$6,$7,$4,$8 }' > \
        !{name}/!{name}.pairs.tmp

    cooler csort -c1 2 -c2 4 \
                 -p1 3 -p2 5 \
                 -p !{task.cpus} \
                 !{name}/!{name}.pairs.tmp \
                 !{chromSizeFile}

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

  publishDir  path: "${params.outputDir}/${name}/matrices/",
              mode: 'copy',
              overwrite: 'true',
              patter: "*.hic"

  input:
  tuple val(name), file(pairs), file(pairsIndex) from resultsPairixJuicer
  file(chromSizeFile) from chromSizeChannelJuicer

  output:
  file("${name}/${name}.hic") into resultsJuicer

  shell:
  juicerGenomes = ['hg18', 'hg19', 'hg38', 'dMel',
                   'mm9', 'mm10', 'anasPlat1', 'bTaurus3',
                   'canFam3', 'equCab2', 'galGal4', 'Pf3D7',
                   'sacCer3', 'sCerS288c', 'susScr3', 'TAIR10']
  genome = juicerGenomes.contains(params.genome) ? params.genome : chromSizeFile
  juicerPath = "${NXF_ASSETS}/t-neumann/hicer-nf/bin
  '''
  mkdir -p !{name}

  java -Xmx!{task.memory} -jar !{juicerPath}/juicer_tools_1.22.01.jar pre \
       -r !{resolutions} \
       -k KR,GW_KR \
       !{pairs} \
       !{name}/!{name}.hic \
       !{genome}
  '''
}

process matrixBuilder {

    tag { name }

    input:
    tuple val(name), file(pairs), file(pairsIndex) from resultsPairixCooler
    file(chromSizeFile) from chromSizeChannelCooler

    output:
    tuple val(name), file("${name}/${name}_1kb.cool") into resultsMatrixBuilder

    shell:
    '''
    mkdir -p !{name}

    cooler cload pairix --assembly !{params.genome} \
                        -p !{task.cpus} \
                        !{chromSizeFile}:1000 \
                        !{pairs} \
                        !{name}/!{name}_1kb.cool
    '''

}

process zoomifyMatrix {

    tag { name }

    input:
    tuple val(name), file(basematrix) from resultsMatrixBuilder

    output:
    tuple val(name), file("${name}/${name}.mcool") into resultsZoomifyMatrix

    shell:
    '''
    mkdir -p !{name}

    cooler zoomify -p !{task.cpus} \
                   -r !{resolutions} \
                   -o !{name}/!{name}.mcool \
                   !{basematrix}
    '''
}

process matrixNormalizer {

    publishDir  path: "${params.outputDir}/${name}/matrices/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "*.mcool"

    tag { name }

    input:
    tuple val(name), file(mcool) from resultsZoomifyMatrix

    output:
    tuple val(name), file("${mcool}") into resultsMatrixNormalizer

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
