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
--samples '[path to samples file]'
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

* [Bowtie 2](https://www.ncbi.nlm.nih.gov/pubmed/22388286/)
  > Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923. PubMed PMID: 22388286; PubMed Central PMCID: PMC3322381.

* [HiCUP](https://www.ncbi.nlm.nih.gov/pubmed/26835000)
  > Wingett SW, Ewels P, Furlan-Magaril M, Nagano T, Schoenfelder S, Fraser P, Simon Andrews S. HiCUP: pipeline for mapping and processing Hi-C data. F1000Research. 2015 4:1310. doi: 10.12688/f1000research.7334.1

* [HiCExplorer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5768762/)
  > Wingett SW, Ewels P, Furlan-Magaril M, Nagano T, Schoenfelder S, Fraser P, Simon Andrews S. High-resolution TADs reveal DNA sequences underlying genome organization in flies. Nature Communications. 2018 9:189. doi: 10.1038/s41467-017-02525-w. PubMed PMID: 5768762.

* [SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943/)
  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

* [MultiQC](https://www.ncbi.nlm.nih.gov/pubmed/27312411/)
  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

### Python packages

* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

* [pyranges](https://github.com/biocore-ntnu/pyranges)
  > Stovner EB, Sætrom P; PyRanges: efficient comparison of genomic intervals in Python. Bioinformatics. 2020 Feb 1;36(3):918-919. doi: 10.1093/bioinformatics/btz615. PMID: 31373614

* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
  
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37
  
* [HiCMatrix](https://github.com/deeptools/HiCMatrix)

### Software packaging/containerisation tools

* [Bioconda](https://www.ncbi.nlm.nih.gov/pubmed/29967506/)
  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

* [Anaconda](https://anaconda.com)
  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Singularity](https://www.ncbi.nlm.nih.gov/pubmed/28494014/)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

* [Docker](https://www.docker.com/)
