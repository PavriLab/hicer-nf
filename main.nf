#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2019 Tobias Neumann, Daniel Malzl
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

        --outputDir      Directory name to save results to. (Defaults to
                         'results')

        References:
        --genome         Name of reference (hg38, mm10)
        --fasta          Alternatively, path to genome fasta file which will be digested

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

if (params.help) {
    helpMessage()
    exit 0
}

params.hicdigest = params.genome ? params.genomes[ params.genome ].hicdigest ?: false : false
params.bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false

log.info ""
log.info " parameters "
log.info " ======================"
log.info " samples list             : ${params.samples}"
log.info " Genome                   : ${params.genome}"
log.info " Fasta                    : ${params.fasta}"
log.info " Bowtie2 index            : ${params.bowtie2}"
log.info " HiC Digest               : ${params.hicdigest}"
log.info " output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

if (params.bowtie2) {
  bowtie2Index = Channel
      .fromPath(params.star_index, checkIfExists: true)
      .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2}" }
} else if (params.fasta) {
  fastaForBowtie2 = Channel.fromPath(params.fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
} else {
    exit 1, "No reference genome files specified!"
}

if (params.hicdigest) {
  hicdigestIndex = Channel
      .fromPath(params.star_index, checkIfExists: true)
      .ifEmpty { exit 1, "HiCdigest not found: ${params.hicdigest}" }
} else if (params.fasta) {
  fastaForHicdigest = Channel.fromPath(params.fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
} else {
    exit 1, "No reference genome files specified!"
}

if (!params.hicdigest && params.fasta) {
      process makeHicDigest {
          tag "$fasta"

          input:
          file fasta from fastaForHicdigest

          output:
          file ("MboI_restriction_sites.bed") into hicdigestIndex

          shell:
          """
          hicup_digester --genome genome --re1 ^GATC,MboI !{fasta}

          findRestSite --fasta !{fasta} --searchPattern GATC -o MboI_restriction_sites.bed
          """
      }
}

if (!params.bowtie2 && params.fasta) {
      process buildBowtie2Index {
          tag "$fasta"

          input:
          file fasta from fastaForBowtie2

          output:
          file "bowtie2" into bowtie2Index

          shell:
          """
          mkdir -p bowtie2Index

          bowtie2-build -f !{fasta} bowtie2Index/genome --threads !{task.cpus}
          """
      }
}

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .set { samples }

process trim {

    tag { parameters.name }

    input:
    val(parameters) from samples

    output:
    file "*_fastqc.{zip,html}" into fastqcResults
    file "*trimming_report.txt" into trimgaloreResults
    set val("${parameters.name}"), file('*_trimmed_val_1.fq.gz'), file('*_trimmed_val_2.fq.gz') into resultsTrimming

    shell:

    """
    trim_galore --paired \
    --quality 20 \
    --fastqc \
    --illumina \
    --gzip \
    --basename !{parameters.name}_trimmed_1.fq \
    --cores !{task.cpus} \
    !{parameters.read1} \
    !{parameters.read2}
    """
}

process trim {

    tag { parameters.name }

    input:
    val(parameters) from samples

    output:
    file "*_fastqc.{zip,html}" into fastqcResults
    file "*trimming_report.txt" into trimgaloreResults
    set val("${parameters.name}"), file('*_trimmed_val_1.fq.gz'), file('*_trimmed_val_2.fq.gz') into resultsTrimming

    shell:

    """
    trim_galore --paired \
    --quality 20 \
    --fastqc \
    --illumina \
    --gzip \
    --basename !{parameters.name}_trimmed_1.fq \
    --cores !{task.cpus} \
    !{parameters.read1} \
    !{parameters.read2}
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
