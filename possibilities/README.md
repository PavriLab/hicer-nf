# Possibilities for downstream analysis
The jupyter notebook in this section shows a possible downstream analysis of high sequencing depth data in which we infer subcompartments from Hi-C data only as described in [Rao et al. 2014](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4).

## Performing clustering of bins using a Gaussian HMM
In order to do so we first need to perform clustering as described using a Gaussian HMM. This can be done with the `performClustering.py` utility which takes a KR-normalized genome-wide matrix, transforms it into a clustering matrix as described and trains a Gaussian HMM to infer a segmentation of the genome for different numbers of clusters.
A usual command to invoke the clustering would be as follows:
```
performClustering.py -m CH12_HiC_200kb_KR.h5 --mink 1 --maxk 15 -r 0.3 -o CH12_HiC_200kb_KR_cluster.npz --imputerows 7689:7905 --imputecols 9592:10116 -pd plots
```

where `--mink` and `--maxk` denote the minimum and maximum number of clusters for which the AIC/BIC should be computed, `-r` specifies the fraction of 0 entries a given row is allowed to have (otherwise it is removed from the matrix before clustering; this value is usually set such that at most 5 - 10% of rows are removed). `--imputerows` and `--imputecols` specify a submatrix in the genome-wide matrix for which the values should be exchanged by imputed values (this is mainly necessary if there are small regions that show unusually strong interchromosomal interaction between an even and an odd chromosome, which would otherwise interfer with the training of the HMM). `-pd` specifies the directory to which the generated clustering visualizations are written.

## Genome-wide interchromosomal contact normalization
One of the downstream analysis steps of the subcompartment inference is the computation of correspondence between clusters on the even chromsosomes and clusters on the odd chromosomes. This is done by computing the enrichment of contacts between bins of a cluster on the even chromosomes and a cluster on the odd chromosomes. In order to remove biases introduced by intrachromosomal contacts (which are naturally more frequent than interchromosomal contacts) we construct a special KR-normalized matrix in which we remove intrachromosomal contacts by setting them to 0 before KR-normalization. To do so we can use the `interchromosomalKRnorm.py` script which is invoked as follows:
```
interchromosomalKRnorm.py h5 -m CH12_HiC_200kb_raw.h5 -o correctedHiC/CH12_HiC_200kb_interKR.npz
```

## Python Packages
