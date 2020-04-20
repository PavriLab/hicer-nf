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

        --resolution     Resolution of the matrix in kb (default: 200)

        --outputDir      Directory name to save results to. (Defaults to
                         'results')

        References:
        --genome         Name of reference (hg38, mm10)
        --fasta          Alternatively, path to genome fasta file which will be digested
        --bowtie2        Optional: Path to bowtie2 index
        --bed12          Genes
        --hicdigest      Restriction site digest index

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
params.hicRestriction = params.genome ? params.genomes[ params.genome ].hicRestriction ?: false : false
params.bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false

log.info ""
log.info " parameters "
log.info " ======================"
log.info " samples list             : ${params.samples}"
log.info " Genome                   : ${params.genome}"
log.info " Fasta                    : ${params.fasta}"
log.info " Bowtie2 index            : ${params.bowtie2}"
log.info " Genes bed                : ${params.bed12}"
log.info " HiC Digest               : ${params.hicdigest}"
log.info " HiC Restriction          : ${params.hicRestriction}"
log.info " output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

if (params.bowtie2) {

  lastPath = params.bowtie2.lastIndexOf(File.separator)
  bwt2_dir =  params.bowtie2.substring(0,lastPath+1)
  bwt2_base = params.bowtie2.substring(lastPath+1)

  bowtie2Index = Channel.fromPath( bwt2_dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bowtie2}" }

} else if (params.fasta) {
  lastPath = params.fasta.lastIndexOf(File.separator)
  bwt2_base = params.fasta.substring(lastPath+1)

  fastaForBowtie2 = Channel.fromPath(params.fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
} else {
    exit 1, "No reference genome files specified!"
}

if (params.hicdigest && params.hicRestriction) {
  Channel
      .fromPath(params.hicdigest, checkIfExists: true)
      .ifEmpty { exit 1, "HiCdigest not found: ${params.hicdigest}" }
      .into{hicdigestIndex ; hicdigestIndexForInsertSize}
  hicRestrictionIndex = Channel
      .fromPath(params.hicRestriction, checkIfExists: true)
      .ifEmpty { exit 1, "HiCRestriction not found: ${params.hicRestriction}" }
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
          file ("Digest*.txt") into hicdigestIndex, hicdigestIndexForInsertSize
          file ("MboI_restriction_sites.bed") into hicRestrictionIndex

          shell:
          """
          hicup_digester --genome genome --re1 ^GATC,MboI !{fasta}

          findRestSite --fasta !{fasta} --searchPattern GATC -o MboI_restriction_sites.bed
          """
      }
}

if (!params.bowtie2 && params.fasta) {
      process buildBowtie2Index {
          tag "$bwt2_base"

          input:
          file fasta from fastaForBowtie2

          output:
          file "bowtie2Index" into bowtie2Index

          shell:
          bwt2_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
          """
          mkdir bowtie2Index

	        bowtie2-build ${fasta} bowtie2Index/${bwt2_base} --threads !{task.cpus}
	        """

      }
}

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .into { samplesChannel ; optionalDiscoveryChannel }

process trim {

    tag { parameters.name }

    input:
    val(parameters) from samplesChannel

    output:
    file "*_fastqc.{zip,html}" into fastqcResults
    file "*trimming_report.txt" into trimgaloreResults
    set val("${parameters.name}"), file('*_trimmed_val_1.fq.gz'), file('*_trimmed_val_2.fq.gz') into resultsTrimming

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

     publishDir path: "${params.outputDir}/QC/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "*/*html"

    input:
    file index from bowtie2Index.collect()
    file digest from hicdigestIndex.collect()
    set val(name), file(fastq1), file(fastq2) from resultsTrimming

    output:
    set val(name), file("${name}/*sam") into resultsHicup, resultsHicupForInsertSize
    file("${name}/*html") into htmlHicup
    file("${name}/HiCUP_summary_report*") into multiqcHicup

    shell:

    '''
    mkdir -p !{name}

    hicup \
    --bowtie2 $(which bowtie2) \
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

optionalDiscoveryChannel
    .map { row -> def flag = true ; if (row.insertSizeMin == null || row.insertSizeMin == '' || row.insertSizeMax == null || row.insertSizeMax == null){ flag = false}; [row.name,row,flag] }
    .set { insertSizeChannel }

process insertSize {

    tag { name }

     publishDir path: "${params.outputDir}/${name}/QC/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "*insertSize_histogram.txt"

    input:
    set val(name), val(parameters), val(insertSize), file(sam) from insertSizeChannel.join(resultsHicupForInsertSize)
    file digest from hicdigestIndexForInsertSize.collect()

    output:
    file("*insertSize_histogram.txt") optional true into histogramInsertSize
    set env(MIN), env(MAX) into resultsInsertSize

    shell:
    if (!insertSize) {
      '''
      getInsertSizeInterval.py -i !{sam} -d !{digest} -o !{name}_insertSize_histogram.txt
      MIN=$(grep "#minInsertSize" !{name}_insertSize_histogram.txt | sed -e "s/.*\\s//g")
      MAX=$(grep "#maxInsertSize" !{name}_insertSize_histogram.txt | sed -e "s/.*\\s//g")
      '''
    } else {
      '''
      MIN=!{parameters.insertSizeMin}
      MAX=!{parameters.insertSizeMax}
      '''
    }
}

process bamPreparation {

    tag { name }

    publishDir path: "${params.outputDir}/${name}/bam/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "*bam"

    input:
    set val(name), file(sam) from resultsHicup

    output:
    set val(name), file("*first.bam"), file("*second.bam") into resultsSamToBam

    shell:

    '''
    samtools view -@ !{task.cpus} -b -f 65 !{sam} > !{name}.hicup.first.bam

    samtools view -@ !{task.cpus} -b -f 129 !{sam} > !{name}.hicup.second.bam

    '''
}

process matrixBuilder {

    tag { name }

    input:
    set val(name), file(first), file(second) from resultsSamToBam
    set env(MIN), env(MAX) from resultsInsertSize

    output:
    set val(name), file("*.h5") into resultsMatrixBuilder

    shell:

    '''

    hicBuildMatrix -s !{first} !{second} -o !{name}_base.h5 --skipDuplicationCheck --binSize 1000 --QCfolder hicQC -ga mm9 --minDistance $MIN --maxLibraryInsertSize $MAX --threads !{task.cpus}

    '''
}

process matrixSubsetter {

    tag { name }

    input:
    set val(name), file(matrix) from resultsMatrixBuilder

    output:
    set val(name), file("*canonical.h5") into resultsMatrixSubsetter, mcoolBuilderProcess

    shell:

    '''
    hicAdjustMatrix -m !{matrix} --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY --action keep -o !{name}_!{params.resolution}kb_canonical.h5
    '''
}

process matrixResolutioner {

  publishDir path: "${params.outputDir}/${name}/matrices/",
         mode: 'copy',
         overwrite: 'true',
         pattern: "*kb.h5"

  tag { name }

  input:
  set val(name), file(matrix) from resultsMatrixSubsetter

  output:
  set val(name), file("*kb.h5") into resultsMatrixResolutioner

  shell:

  '''
  hicMergeMatrixBins -m !{matrix} -o !{name}_!{params.resolution}kb.h5 -nb !{params.resolution}

  '''
}

process matrixNormalizer {

    publishDir path: "${params.outputDir}/${name}/matrices/",
             mode: 'copy',
             overwrite: 'true',
             pattern: "*KR.h5"

    tag { name }

    input:
    set val(name), file(matrix) from resultsMatrixResolutioner

    output:
    set val(name), file("*canonical_KR.h5") into resultsMatrixNormalizer

    shell:

    '''
    hicCorrectMatrix correct -m !{matrix} --correctionMethod KR -o !{name}_!{params.resolution}kb_canonical_KR.h5
    '''
}

process matrixEO {

    publishDir path: "${params.outputDir}/${name}/matrices/",
             mode: 'copy',
             overwrite: 'true',
             pattern: "*EO.h5"

    tag { name }

    input:
    set val(name), file(matrix) from resultsMatrixNormalizer

    output:
    set val(name), file("*canonical_EO.h5") into resultsMatrixEO

    shell:

    '''
    hicTransform -m !{matrix} --method obs_exp_lieberman -o !{name}_!{params.resolution}kb_canonical_EO.h5
    '''
}

process compartmentalization {

    publishDir path: "${params.outputDir}/${name}/bigwigs/",
             mode: 'copy',
             overwrite: 'true',
             pattern: "*bw"

    tag { name }

    input:
    set val(name), file(matrix) from resultsMatrixEO

    output:
    set val(name), file("*.bw") into resultsCompartmentalization

    shell:

    '''
    hicAdjustMatrix -m !{matrix} --chromosomes chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY --action keep -o tmp.h5

    generateEigenvectorBigWig.py -m tmp.h5 -g !{params.bed12} -r !{params.resolution.toInteger() * 1000} --chromLengths /groups/zuber/zubarchive/USERS/tobias/mm9/mm9.chr_lengths.txt -o !{name}_!{params.resolution}kb_eigv.bw
    '''
}

process mcoolBuilder {

  tag { name }

  publishDir path: "${params.outputDir}/${name}/mcool/",
           mode: 'copy',
           overwrite: 'true',
           pattern: "*.mcool"

  input:
  set val(name), file(matrix) from mcoolBuilderProcess

  output:
  set val(name), file("*.mcool") into resultsmcoolBuilder

  shell:

  '''
  for i in 5 10 25 50 100 500 1000;
  do
    hicMergeMatrixBins -m !{matrix} -o !{name}_${i}kb.h5 -nb ${i}
    hicCorrectMatrix correct -m !{name}_${i}kb.h5 --correctionMethod KR -o !{name}_${i}kb_KR.h5
  done

  h5toCool.py -m *KR.h5 -o !{name}.mcool --merge

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
