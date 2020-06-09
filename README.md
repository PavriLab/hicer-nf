# hicer-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

**hicer-nf** is a bioinformatics analysis pipeline used for [HiC](https://en.wikipedia.org/wiki/Chromosome_conformation_capture#Hi-C_(all-vs-all)) data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline summary

We loosely follow the steps proposed by [Rao et al, Cell 2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547) for processing HiC datasets.

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment, Filtering and Deduplication ([`HiCUP`](https://www.bioinformatics.babraham.ac.uk/projects/hicup/))
4. Generating an indexed pairs file from the alignments ([`pairix`](https://github.com/4dn-dcic/pairix))
5. Build base resolution contact matrix (i.e. the smallest resolution given; default 5kb) ([`cooler`](https://cooler.readthedocs.io/en/latest/))
6. Aggregate bins to a range of default resolutions including optional custom resolutions ([`cooler`](https://cooler.readthedocs.io/en/latest/))
7. Compute matrix normalization vectors for all aggregated resolutions with the [KR](https://doi.org/10.1093/imanum/drs019) and the [IC](https://www.nature.com/articles/nmeth.2148) algorithm (HiCExplorers [C++ implementation of the KR algorithm](https://github.com/deeptools/Knight-Ruiz-Matrix-balancing-algorithm) and [`cooler`](https://cooler.readthedocs.io/en/latest/) Out-Of-Core balancing)
8. Present QC for raw reads, alignment and filtering [`MultiQC`](http://multiqc.info/)

The generated mcool files are compatible with [cooltools](https://cooltools.readthedocs.io/en/latest/index.html) for downstream analysis and [Higlass](https://github.com/higlass/higlass) for visualization. A range of possible downstream analyses can be found in the possibilities directory of this repo. In addition the pipeline also generates [hic] files with [juicer_tools](https://github.com/aidenlab/juicer/wiki/Download) pre for visualization in [juicebox](https://github.com/aidenlab/Juicebox) and/or loop calling with juicer's HICCUPS algorithm. 

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Start running your own analysis!

```bash
# if you have specified an igenomes directory from or run a profile from any of the institutions 
# listed here https://raw.githubusercontent.com/nf-core/configs/master/nfcore_custom.config
nextflow run t-neumann/hicer-nf --samples samples.txt --genome mm9 --re1 ^GATC,MboI

# else you can specify the most essential things manually
# if HICUP digest and bowtie2 index are not available
nextflow run t-neumann/hicer-nf --samples samples.txt --genome mm9 --re1 ^GATC,MboI --fasta genome.fa --chromSizes chrom.sizes

# if HICUP digest and bowtie2 index are available
nextflow run t-neumann/hicer-nf --samples samples.txt --genome mm9 --re1 ^GATC,MboI --bowtie2Index /path/to/bowtie2Index/genome_base_name --hicupDigest /path/to/hicup/digest/file --chromSizes chrom.sizes
```

These invocations compute cooler and hic files for a default resolution list of 5kb, 10kb, 25kb, 50kb, 100kb, 250kb, 500kb and 1Mb. If you want resolutions that are not listed here you could use the `--resolutions` parameter (see below).

## Main arguments

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. The pipeline uses the [nf-core custom_configs repo](https://github.com/nf-core/configs) to load config files for compute clusters of different institutions. So please check with this or make your own config file if you want to use the igenomes reference database. Note that multiple profiles can be loaded, for example: `-profile docker cbe` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`zuberlab/hicer-nf`](http://hub.docker.com/r/zuberlab/hicer-nf/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`zuberlab/hicer-nf`](http://hub.docker.com/r/zuberlab/hicer-nf/)

#### `--samples`

You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a tab-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--samples '[path to samples file]'
```

```bash

name read1 read2
WT WT_1.fastq.gz WT_2.fastq.gz
KD KD_1.fastq.gz KD_2.fastq.gz
```

| Column      | Description                                                                                                 |
|-------------|-------------------------------------------------------------------------------------------------------------|
| `name` | Name of this sample.
| `read1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `read2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |

#### `--genome`

The genome argument has two properties. Firstly, it is used to document the name of the reference genome of the organism the Hi-C data originates from and secondly it can be used to retrieve prespecified reference data from a local igenomes database, where the pipeline automatically takes the files it requires for processing the Hi-C data (i.e. a bowtie2 index, a genome fasta and a chromSizes file).

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta      = "/path/to/genome/fasta/file" // Used if no bowtie2 index or hicupDigest given
      bowtie2    = "/path/to/bowtie2/index/basename
      chromSizes = "/path/to/chrom/sizes/file
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

#### `--re1`

Hi-C experiments typically involve digestion of the with restriction enzymes after cross-linking. To be able to computationally identify artifacts in the sequence data arising from this process the HICUP pipeline requires the sequence motif the used enzyme cuts including its name. This information is given via the `--re1` parameter in the form specified by [HICUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html) which is ^GATC,MboI, where ^ indicates the cutsite, GATC is the sequence motif the enzyme recognizes and MboI is the name of the enzyme. In cases where a double digest protocol was used one can also specify `--re2` for the second restriction enzyme used (see below).

#### `--fasta`

This parameter is used to specify the genome fasta file and is only required if no igenomes database is available or hicup digestion file or bowtie2 index is not available locally. The file is used for in-silico restriction digestion for HICUP (if the file is not specified manually with `--hicupDigest`, see below) and bowtie2 index computation (if not specified manually)

#### `--chromSizes`

This parameter is used to specify file containing the chromosome names and their size. This information is used by cooler for binning the genome. The file can either be a tab-separated file with two columns (preferred format)

```bash
cat chrom.sizes.tsv

chr1 10000
chr2 30499
```

or an XML file as given by the igenomes database

```bash
cat chrom.sizes.xml

<sequenceSizes genomeName="genome">
        <chromosome fileName="genome.fa" contigName="chr1" totalBases="10000" isCircular="false" md5="be7e6a13cc6b9da7c1da7b7fc32c5506" ploidy="2" knownBases="126847849" />
        <chromosome fileName="genome.fa" contigName="chr2" totalBases="30499" isCircular="false" 
</sequenceSizes>
```

The XML will be converted to TSV in-situ where the pipeline uses the `contigName` and `totalBases` keys for generating the TSV file. Thus at least these fields have to be present in the XML. This option is for compatibility with older igenomes databases.

## Generic arguments

#### `--bowtie2Index`

Full path to an existing bowtie2 index for your reference genome including the base name for the index. This is only necessary if you don't have an igenomes database and have a precomputed index for your fasta file that you want to use. Otherwise, the index will be computed from the fasta file specified with `--fasta`

```bash
--bowtie2 '[directory containing bowtie2 index]/genome'
```

#### `--hicupDigest`

File with digested genome of a given restriction enzyme as produced with `hicup_digester`. This can be used to override in-situ digestion of the fasta and use your own digestion file.

```bash
--hicdigest '[file with digested genome]'
```

#### `--resolutions`

By default, the pipeline computes matrices and normalizations thereof for resolutions 5kb, 10kb, 25kb, 50kb, 100kb, 250kb, 500kb and 1Mb for both cooler and hic files. However, if you want to add additional resolutions which are not present in this list you can use the `--resolutions` parameter by passing either a single value or multiple ones separated by commas in Bp. The specified resolutions will then be added. However, be aware that the the specified resolutions have to be an integer multiple of the smallest resolution in the list i.e. a multiple of 5kb when using the default resolutions or a multiple of any resolution specified that is smaller than 5kb. In this sense any smaller resolution should also be a divisor of all default resolutions (e.g. 500bp, 1kb would all be fine however 1.5kb will not work because 5kb is not an integer multiple of 1.5kb)

```bash
--resolution 20000 / --resolution 20000,75000,200000
```

#### `--outputDir`

Name of the folder to which the output will be saved (default: results)

```bash
--outputDir '[directory name]'
```

## Credits

The pipeline was developed by [Tobias Neumann](mailto:tobias.neumann.at@gmail.com) and [Daniel Malzl](mailto:daniel.malzl@gmx.at) for use at the [IMP](https://www.imp.ac.at/), Vienna.

The [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [nf-core/chipseq](https://github.com/nf-core/chipseq) pipelines developed by Phil Ewels were initially used as a template for this pipeline. Many thanks to Phil for all of his help and advice, and the team at SciLifeLab.

Many thanks to others who have helped out along the way too, including (but not limited to): [@apeltzer](https://github.com/apeltzer), [@micans](https://github.com/micans), [@pditommaso](https://github.com/pditommaso).

## Citations

### Pipeline tools

* [Nextflow](https://www.ncbi.nlm.nih.gov/pubmed/28398311/)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

* [Bowtie 2](https://www.ncbi.nlm.nih.gov/pubmed/22388286/)
  > Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923. PubMed PMID: 22388286; PubMed Central PMCID: PMC3322381.
  
* [Juicer](https://www.cell.com/cell-systems/fulltext/S2405-4712%2816%2930219-8)
  > Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell Systems 3(1), 2016. doi: 10.1016/j.cels.2016.07.002

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

* [pairix](https://github.com/4dn-dcic/pairix)

### Python packages

* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
  
* [cooler](https://cooler.readthedocs.io/en/latest/)
  > Nezar Abdennur, Leonid A Mirny, Cooler: scalable storage for Hi-C data and other genomically labeled arrays, Bioinformatics, Volume 36, Issue 1, 1 January 2020, Pages 311–316. doi: 10.1093/bioinformatics/btz540
  
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37
  
* [scipy](https://www.scipy.org/)
  > Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). doi: 10.1038/s41592-019-0686-2
  
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.
  
* [h5py](https://github.com/h5py/h5py)

### Software packaging/containerisation tools

* [Bioconda](https://www.ncbi.nlm.nih.gov/pubmed/29967506/)
  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

* [Anaconda](https://anaconda.com)
  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Singularity](https://www.ncbi.nlm.nih.gov/pubmed/28494014/)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

* [Docker](https://www.docker.com/)
