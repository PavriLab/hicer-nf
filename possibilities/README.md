# Possibilities for downstream analysis
The jupyter notebook in this section shows a possible downstream analysis of high sequencing depth data in which we infer subcompartments from Hi-C data only as described in [Rao et al. 2014](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4).

## Performing clustering of bins using a Gaussian HMM
In order to do so we first need to perform clustering as described using a Gaussian HMM. This can be done with the `performClustering.py` utility which takes a KR-normalized genome-wide matrix, transforms it into a clustering matrix as described and trains a Gaussian HMM to infer a segmentation of the genome for different numbers of clusters.
A usual command to invoke the clustering would be as follows:
```
performClustering.py -m CH12_HiC_200kb_KR.h5 --mink 1 --maxk 15 -r 0.3 -o CH12_HiC_200kb_KR_cluster.npz --imputerows 7689:7905 --imputecols 9592:10116 -pd plots
```

where `--mink` and `--maxk` denote the minimum and maximum number of clusters for which the clustering should be computed, `-r` specifies the fraction of 0 entries a given row is allowed to have (otherwise it is removed from the matrix before clustering; this value is usually set such that at most 5 - 10% of rows are removed). `--imputerows` and `--imputecols` specify a submatrix in the genome-wide matrix for which the values should be exchanged by imputed values (this is mainly necessary if there are small regions that show unusually strong interchromosomal interaction between an even and an odd chromosome, which would otherwise interfer with the training of the HMM). `-pd` specifies the directory to which the generated clustering visualizations are written.
The script performs model training and inference of the most probable sequence of clusterassignments over all rows in the cluster matrix and generates multiple result files:

1.  `*_cluster.npz`
This file is a compressed collection of numpy arrays and contains several components as follows:
```
/
 |── evenk${mink}; numpy.array holding bin clusterassignments of even chromosomes for HMM with k components
 :        :
 |── evenk${maxk}
 |
 |── oddk${mink}; numpy.array holding bin clusterassignments of odd chromosomes for HMM with k components
 :        :
 |── oddk${maxk}
 |
 |── evenremrows; numpy.array holding indices of rows that were removed from the matrix during generation of the even chromosome clustermatrix
 |
 |── evenremcols; numpy.array holding indices of columns that were removed from the matrix during generation of the even chromosome clustering matrix
 |
 |── oddremrows; numpy.array holding indices of rows that were removed from the matrix during generation of the odd chromosome clustermatrix
 |
 |── oddremcols; numpy.array holding indices of columns that were removed from the matrix during generation of the odd chromosome clustermatrix
 |
 |── evenAIC; numpy.array holding AIC values for all models from mink to maxk for even chromosomes
 |
 |── evenBIC; numpy.array holding BIC values for all models from mink to maxk for even chromosomes
 |
 |── oddAIC; numpy.array holding AIC values for all models from mink to maxk for odd chromosomes
 |
 └── oddBIC; numpy.array holding BIC values for all models from mink to maxk for odd chromosomes
    
```
The contained data can easily be retrieved using `numpy.load`

2. `*_informationcriterion.pdf`
Visualization of AIC and BIC for even and odd clustering matrices (is saved to the directory specified with `-pd`)

3. `*_[even|odd]_k${k}.png`
Visualization of the different clustering results for the clustering matrices where `${k}` is between `--mink` and `--maxk` (is saved to the directory specified with `-pd`)

## Genome-wide interchromosomal contact normalization
One of the downstream analysis steps of the subcompartment inference is the computation of correspondence between clusters on the even chromsosomes and clusters on the odd chromosomes. This is done by computing the enrichment of contacts between bins of a cluster on the even chromosomes and a cluster on the odd chromosomes. In order to remove biases introduced by intrachromosomal contacts (which are naturally more frequent than interchromosomal contacts) we construct a special KR-normalized matrix in which we remove intrachromosomal contacts by setting them to 0 before KR-normalization. To do so we can use the `interchromosomalKRnorm.py` script which is invoked as follows:
```
interchromosomalKRnorm.py -m CH12_HiC_200kb_raw.h5 -o correctedHiC/CH12_HiC_200kb_interKR.npz
```

The generated file `*_interKR.npz` is similar to the clustering matrices where rows are bins on odd chromosomes and columns are bins on even chromosomes.

## Python Packages
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37

* [scipy](https://www.scipy.org/)
  > Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). doi: 10.1038/s41592-019-0686-2
  
* [pytables](https://www.pytables.org/)

* [hmmlearn](https://github.com/hmmlearn/hmmlearn)

* [matplotlib](https://matplotlib.org/)
  > John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007). doi: 10.1109/MCSE.2007.55

* [krbalancing](https://github.com/deeptools/Knight-Ruiz-Matrix-balancing-algorithm)
