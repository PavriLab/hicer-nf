# Possibilities for downstream analysis
The jupyter notebooks in this section show some possible downstream analysis. Specifically, generic analyses and special ones like of high sequencing depth data in which we infer subcompartments from Hi-C data only as described in [Rao et al. 2014](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4).

## Plotting contact matrices
The `plotHiCmat.py` script is a lightweight utility to visualize your generated contact matrices. An example command would be:
```
plotHiCmat.py -m CH12_HiC_200kb_KR.h5 --vMax 50 -o CH12_HiC_200kb_plot.jpg
```

The `--vMax` specifies the maximum value for the colormap used to plot the matrix (`--vMin` can be used to set the minimum value). In addition to this one can also specify the chromosomes that should be included in the plot with the `--chromosomes` argument, which takes a space-separated list of chromosomes. The used colormap can be changed with the `--colorMap` argument (this has to be either `redmap` (default) or one of the named colormaps offered by `matplotlib`). 

## Performing clustering of bins using a Gaussian HMM
In order to perform subcompartment inference we first have to convert the a resolution of the mcooler file to h5 format using [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/)(note that 200kb is not part of the default resolutions but you can always compute it post-pipeline using `cooler coarsen`)

```
hicConvertFormat -m mcoolfile.mcool::resolutions/200000 \
                 --inputFormat cool \
                 --outputFormat h5 \
                 --correction_name KR \
                 -o 200kb.h5
```

We can then perform clustering as described using a Gaussian HMM. This can be done with the `performClustering.py` utility which takes a KR-normalized genome-wide matrix, transforms it into a clustering matrix as described and trains a Gaussian HMM to infer a segmentation of the genome for different numbers of clusters. A usual command to invoke the clustering would be as follows:
```bash
performClustering.py -m CH12_HiC_200kb_KR.h5 \
                     --mink 1 \
                     --maxk 15 \
                     -r 0.3 \
                     -p CH12_HiC_200kb_KR_cluster \
                     --imputerows 7689:7905 \
                     --imputecols 9592:10116 \
                     -pd plots
```

where `--mink` and `--maxk` denote the minimum and maximum number of clusters for which the clustering should be computed, `-r` specifies the fraction of 0 entries a given row is allowed to have (otherwise it is removed from the matrix before clustering; this value is usually set such that at most 5 - 10% of rows are removed). `--imputerows` and `--imputecols` specify a submatrix in the genome-wide matrix for which the values should be exchanged by imputed values (this is mainly necessary if there are small regions that show unusually strong interchromosomal interaction between an even and an odd chromosome, which would otherwise interfer with the training of the HMM). `-pd` specifies the directory to which the generated clustering visualizations are written. The script performs model training and inference of the most probable sequence of clusterassignments over all rows in the cluster matrix and generates multiple result files:

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

4. `*.hmm*.pkl`
These files contain the trained model for each number of states chosen, dumped with `pickle` and can be loaded with
```
import pickle
model = pickle.load('*.hmm*.pkl')
```

## Genome-wide interchromosomal contact normalization
One of the downstream analysis steps of the subcompartment inference is the computation of correspondence between clusters on the even chromsosomes and clusters on the odd chromosomes. This is done by computing the enrichment of contacts between bins of a cluster on the even chromosomes and a cluster on the odd chromosomes. In order to remove biases introduced by intrachromosomal contacts (which are naturally more frequent than interchromosomal contacts) we construct a special KR-normalized matrix in which we remove intrachromosomal contacts by setting them to 0 before KR-normalization. To do so we can use the `interchromosomalKRnorm.py` script which is invoked as follows:
```
interchromosomalKRnorm.py -m CH12_HiC_200kb_raw.h5 -o correctedHiC/CH12_HiC_200kb_interKR.npz
```

The generated file `*_interKR.npz` is similar to the clustering matrices where rows are bins on odd chromosomes and columns are bins on even chromosomes.

## Further downstream analysis of the clustering results
The notebook `subcompartment.ipynb` contains further downstream analyses of the obtained clustering to generate a genome-wide annotation of subcompartments.

## Generic downstream analyses
The notebook `hicanalysis.ipynb` contains generic downstream analyses schemes using the [cooltools](https://cooltools.readthedocs.io/en/latest/index.html) package. Here you find things like making aggregates of loops and TADs, computing insulation scores and making saddle plots. Eigenvectors can typically be optained with [HOMER](http://homer.ucsd.edu/homer/interactions/)(as in our case) or other packages like [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/). However, be aware that cooltools is based on the [cooler](https://cooler.readthedocs.io/en/latest/index.html) and by the time of this writing the [`cooler.matrix`](https://cooler.readthedocs.io/en/latest/api.html#cooler.Cooler.matrix) function assumes that KR balancing weights are divisive as anticipated from the juicer_tools implementation and [`hic2cool`](https://github.com/4dn-dcic/hic2cool) does not change this anymore. In contrast, the HiCExplorer implementation of the KR algorithm yields multiplicative weights thus if you are using the KR with cooltools make sure that `divisive_weights = False` if `balance = True`) like

```python
clr = cooler.Cooler('sample.mcool::resolutions/100000')
clr.matrix(balance = True, divisive_weights = False)
```

otherwise cooler will balance the matrix by applying the KR weights in a divisive manner and the results will be very large and incorrect. Thus, you might have a look at the cooltools code and implement the functions accordingly for the tasks you want to achieve. The `hicanalysis.ipynb` jupyter notebook contains a guide in code (and a bit of comments) to get a typical downstream analysis started with the results of the pipeline.

## Python Packages
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37

* [scipy](https://www.scipy.org/)
  > Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). doi: 10.1038/s41592-019-0686-2

* [matplotlib](https://matplotlib.org/)
  > John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007). doi: 10.1109/MCSE.2007.55
  
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
  
* [pytables](https://www.pytables.org/)

* [hmmlearn](https://github.com/hmmlearn/hmmlearn)

* [krbalancing](https://github.com/deeptools/Knight-Ruiz-Matrix-balancing-algorithm)

* [cooltools](https://cooltools.readthedocs.io/en/latest/index.html)

* [pypairix](https://pypi.org/project/pypairix/)
