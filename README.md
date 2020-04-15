# hicer-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

**hicer-nf** is a bioinformatics analysis pipeline used for [HiC](https://en.wikipedia.org/wiki/Chromosome_conformation_capture#Hi-C_(all-vs-all)) data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline summary

We loosely follow the steps proposed by [Rao et al, Cell 2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547) for processing HiC datasets.

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment ([`HiCUP`](https://www.bioinformatics.babraham.ac.uk/projects/hicup/))
4. Insert size calculation
5. Build 1kb resolution contact matrix ([`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/))
6. Merge bins to desired resolution and remove non-canonical chromosomes ([`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/))
7. [KR normalize](https://doi.org/10.1093/imanum/drs019) contact matrix ([`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/))
8. Calculate observed / expected matrix ([`HiCExplorer`](https://hicexplorer.readthedocs.io/en/latest/))
9. Create Eigenvector bigwig tracks for A/B compartmentalization ([`pyBigWig`](https://github.com/deeptools/pyBigWig))
10. Create `mcooler` files for use in [HiGlass](http://higlass.io/) ([`cooler`](https://mirnylab.github.io/cooler/))
11. Present QC for raw reads, alignment and filtering [`MultiQC`](http://multiqc.info/)

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Start running your own analysis!

```bash
nextflow run t-neumann/hicer-nf --design design.txt --genome mm9 --singleEnd
```

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`zuberlab/hicer-nf`](http://hub.docker.com/r/zuberlab/hicer-nf/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`zuberlab/repliseq-nf`](http://hub.docker.com/r/zuberlab/hicer-nf/)

### `--samples`

You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--samples '[path to design file]'
```

```bash

name,read1,read2
WT,WT_1.fastq.gz,WT_2.fastq.gz
KD,KD_1.fastq.gz,KD_2.fastq.gz
```

| Column      | Description                                                                                                 |
|-------------|-------------------------------------------------------------------------------------------------------------|
| `name` | Name of this sample.
| `read1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `read2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |

## Generic arguments

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no bwa index given
      bwa     = '<path to the bwa index file>'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--fasta`

Full path to fasta file containing reference genome (*mandatory* if `--genome` is not specified). If you don't have a BWA index available this will be generated for you automatically. Combine with `--save_reference` to save BWA index for future runs.

```bash
--fasta '[path to FASTA reference]'
```

### `--bowtie2`

Full path to an existing Bowtie index for your reference genome including the base name for the index.

```bash
--bowtie2 '[directory containing bowtie2 index]/genome.fa'
```

### `--bed12`

Full path to an existing gene bed file for A/B compartmentalization.

```bash
--bed12 '[bed file with gene coordinates]'
```

### `--hicdigest`

File with digested genome of a given restriction enzyme as produced with `hicup_digester`.

```bash
--hicdigest '[file with digested genome]'
```

### `--resolution`

Resolution of the matrix for KR-normalization, expected/observed calculation and A/B compartment track calculation.

```bash
--resolution '[resolution in kb]'
```

### `--outputDir`

Name of the folder to which the output will be saved (default: results)

```bash
--outputDir '[directory name]'
```

## Credits

The pipeline was developed by [Tobias Neumann](mailto:tobias.neumann.at@gmail.com) for use at the [IMP](https://www.imp.ac.at/), Vienna.

The [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [nf-core/chipseq](https://github.com/nf-core/chipseq) pipelines developed by Phil Ewels were initially used as a template for this pipeline. Many thanks to Phil for all of his help and advice, and the team at SciLifeLab.

Many thanks to others who have helped out along the way too, including (but not limited to): [@apeltzer](https://github.com/apeltzer), [@micans](https://github.com/micans), [@pditommaso](https://github.com/pditommaso).

## Citations

### Pipeline tools

* [Nextflow](https://www.ncbi.nlm.nih.gov/pubmed/28398311/)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

* [BWA](https://www.ncbi.nlm.nih.gov/pubmed/19451168/)
  > Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009 Jul 15;25(14):1754-60. doi: 10.1093/bioinformatics/btp324. Epub 2009 May 18. PubMed PMID: 19451168; PubMed Central PMCID: PMC2705234.

* [BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/20110278/)
  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PubMed PMID: 20110278; PubMed Central PMCID: PMC2832824.

* [SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943/)
  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

* [BamTools](https://www.ncbi.nlm.nih.gov/pubmed/21493652/)
  > Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011 Jun 15;27(12):1691-2. doi: 10.1093/bioinformatics/btr174. Epub 2011 Apr 14. PubMed PMID: 21493652; PubMed Central PMCID: PMC3106182.

* [UCSC tools](https://www.ncbi.nlm.nih.gov/pubmed/20639541/)
  > Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010 Sep 1;26(17):2204-7. doi: 10.1093/bioinformatics/btq351. Epub 2010 Jul 17. PubMed PMID: 20639541; PubMed Central PMCID: PMC2922891.

* [deepTools](https://www.ncbi.nlm.nih.gov/pubmed/27079975/)
  > Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257. Epub 2016 Apr 13. PubMed PMID: 27079975; PubMed Central PMCID: PMC4987876.

* [MultiQC](https://www.ncbi.nlm.nih.gov/pubmed/27312411/)
  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

* [picard-tools](http://broadinstitute.github.io/picard)

### R packages

* [R](https://www.R-project.org/)
  > R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

* [preprocessCore](https://cran.r-project.org/web/packages/getopt/index.html)
  > Bolstad B (2019). preprocessCore: A collection of pre-processing functions.

* [getopt](https://cran.r-project.org/web/packages/getopt/index.html)
  > Trevor L Davis (2010). getopt: C-Like 'getopt' Behavior.

### Software packaging/containerisation tools

* [Bioconda](https://www.ncbi.nlm.nih.gov/pubmed/29967506/)
  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

* [Anaconda](https://anaconda.com)
  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Singularity](https://www.ncbi.nlm.nih.gov/pubmed/28494014/)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

* [Docker](https://www.docker.com/)
