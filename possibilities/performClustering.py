#!/usr/bin/env python

import argparse as ap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from hmmlearn.hmm import GaussianHMM
import scipy.stats as scistats
import logging
import pickle
import os, ntpath
import tables
import cooler
from scipy.sparse import csr_matrix, triu, lil_matrix

def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s

    if isinstance(s, bytes):  # or isinstance(s, np.bytes_):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')

    if isinstance(s, list):
        return [toString(x) for x in s]

    if isinstance(s, np.ndarray):
        return s.astype(str)

    return s


def loadH5(filename, includechroms=None, csr=True, returnintervals = False, dtype = int):
    '''
    loadH5(filename, includechroms=None, csr=True, returnintervals = False)

    loads an *.h5 hic matrix as created by hicexplorer

    :param filename:        name of the *.h5 file containing the matrix
    :param includechroms:   list of chromosomes to include in the returned objects
                            if not given all chromosomes in the *.h5 file are included
    :param csr:             if True returns a csr_matrix object else a full numpy.array
    :param returnintervals: if True also returns the intervals read

    :return:                csr_matrix containing the data in the matrix
    '''
    with tables.open_file(filename) as f:
        parts = {}
        try:
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
        except Exception:
            logging.info('No h5 file. Please check parameters concerning the file type!')
            exit(1)

        matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                            shape=parts['shape'], dtype=dtype)

        intervals = {}
        for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
            if toString(interval_part) == toString('chr_list'):
                chrom_list = getattr(f.root.intervals, interval_part).read()
                intervals[interval_part] = toString(chrom_list)
            else:
                intervals[interval_part] = getattr(f.root.intervals, interval_part).read()

        cut_intervals = list(
            zip(intervals['chr_list'], intervals['start_list'], intervals['end_list'], intervals['extra_list']))

        assert len(cut_intervals) == matrix.shape[0], \
            "Error loading matrix. Length of bin intervals ({}) is different than the " \
            "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

        # compute index array and chromosome list
        inds, chr_list, chroms = [], [], set()
        for i, (chr, start, end, extra) in enumerate(cut_intervals):
            if chr not in chroms:
                chroms.add(chr)
                inds.append(i)
                chr_list.append(chr)

        # if includechroms is given we filter the output for the chromosomes listed
        # and recompute indices of chromosome boundaries in the resulting matrix
        if includechroms:
            includechroms = set(includechroms)
            filterinds, filterchrs = [], []
            for i, chr in zip(range(len(inds)), chr_list):
                if chr in includechroms:
                    filterinds.append([inds[i], inds[i + 1] if i + 1 != len(inds) else matrix.shape[0]])
                    filterchrs.append(chr)

            matrixinds = np.zeros(shape=matrix.shape[0], dtype=bool)
            ncuts, tmpe = [], 0
            for s, e in filterinds:
                matrixinds[s: e] = True

                if s == tmpe:
                    ncuts.append(s)
                    tmpe = e

                else:
                    ncuts.append(tmpe)
                    tmpe = e - s + tmpe


            matrix = matrix[matrixinds, :][:, matrixinds]

            inds = ncuts

            chr_list = filterchrs

    if not csr:
        x = matrix.toarray()
        xi, yi = np.triu_indices(x.shape[0], k=1)
        x[yi, xi] = x[xi, yi]
        matrix = x

    if returnintervals:
        return matrix, np.array(inds), np.array(chr_list), intervals

    else:
        return matrix, np.array(inds), np.array(chr_list)


def loadCooler(cooleruri, applyNorm = False, norm = 'weight', includeChroms = None, nans_to_zero = False):
    '''
    loads a cooler into a csr matrix
    taken from HiCMatrix cool.py see also
    https://github.com/deeptools/HiCMatrix/blob/master/hicmatrix/lib/cool.py

    :param cooleruri:       uri to a given cooler
    :param applyNorm:       if True then the 'norm' is applied to the datapoints in the matrix
    :param norm:            normalization weights to apply if applyNorm is set True
    :param includeChroms:   list of chromosomes to load, if given only the specified chromosomes will be loaded from the cooler

    :return:            data in cooler as scipy.sparse.csr_matrix
    '''
    cooler_file = cooler.Cooler(cooleruri)
    matrix = cooler_file.matrix(balance = norm if applyNorm else False)[:]

    chroms = cooler_file.chromnames
    inds = set()
    for chrom in chroms:
        for binidx in cooler_file.extent(chrom):
            inds.add(binidx)

    inds = sorted(list(inds))

    if includeChroms:
        includechroms = set(includeChroms)
        filterinds, filterchroms = [], []
        for i, chr in zip(range(len(inds)), chroms):
            if chr in includechroms:
                filterinds.append([inds[i], inds[i + 1] if i + 1 != len(inds) else matrix.shape[0]])
                filterchroms.append(chr)

        matrixinds = np.zeros(shape=matrix.shape[0], dtype=bool)
        ncuts, tmpe = [], 0
        for s, e in filterinds:
            matrixinds[s: e] = True

            if s == tmpe:
                ncuts.append(s)
                tmpe = e

            else:
                ncuts.append(tmpe)
                tmpe = e - s + tmpe

        matrix = matrix[matrixinds, :][:, matrixinds]

        inds = ncuts

        chroms = filterchroms

        if nans_to_zero:
            matrix[np.isnan(matrix)] = 0

    return matrix, np.array(inds), np.array(chroms)


def constructClusterContactMatrix(gwmat, chrlist, indarr, excluderows = None, excludecols = None,
                                  imputerows = None, imputecols = None, removelim = 0.3, withX = False,
                                  even = False, transform = True):
    '''
    constructClusterContactMatrix(gwmat, chrlist, indarr, removelim = 0.3, excluderows = None, excludecols = None,
                                  imputerows = None, imputecols = None, withX = False, even = False, transform = True)

    given a normalized, genomewide contact matrix (can be constructed with
    ConstructGenomeWideContactMatrix) constructs a matrix C suitable for performing
    clustering as described in Rao et al. 2014. In particular C is constructed such
    that C_i,j contains the normalized interaction between odd chromosome i and even
    chromosome j. Rows and columns with a number zeros or undefined entries larger than
    removelim of the row/col are removed. Note that the bins to be removed are computed
    sequentially first rows then columns. If even is True the matrix is transposed
    prior to row/col removal to keep removal for odd and even chromosomes consistent
    rows and columns given by excluderows/cols are excluded from the analysis. However,
    you can also pass a list of rows and columns using the imputerows/cols to impute specific
    values specified by row/col with a random value drawn from the rows distribution before
    z-score transformation

    :param gwmat:       genomewide normalized contact matrix
    :param chrlist:     sorted list of chromosomes in gwmat
                        (ascending chr1, chr2, ..., chr10, chr11, ..., chrX
    :param indarr:      array containing the indices of the single matrices in gwmat
                        see ConstructGenomeWideContactMatrix for more details
    :param excluderows: list of row indices corresponding to indarr that should be excluded
    :param excludecols: list of column indices corresponding to indarr that should be excluded
                        from the clustering matrix
    :param imputerows:  list of rows for which an imputation should be performed
    :param imputecols:  list of cols for which an imputation should be performed
    :param removelim:   limit of fraction of undefined or zero entries in row/col
                        rows/cols with sum(0 | NaN)/len(row | col) > 0.3 are removed
    :param withX:       if True chromosome X is included in the even chromosomes
    :param even:        if True Cij is transposed prior to row/col removal
    :param transform:   if True logarithm and zscore transformation are applied to Cij

    :return:            contact subset matrix where rows are only composed of odd chromosomes
                        and columns of even chromosomes (including X if withx = True, or vice versa
                        if even = True) rowindices that were removed, column indices that were removed
    '''

    # building boolean index arrays for row and columns
    colindex = np.zeros(shape = gwmat.shape[1], dtype = bool)
    rowindex = np.zeros(shape = gwmat.shape[0], dtype = bool)

    if excluderows or imputerows:
        processrowchroms = set()
        for indlist, indtype in zip([excluderows, imputerows], ['exclude', 'impute']):
            if not indlist:
                continue

            else:
                rowcounts, rowbins = np.histogram(indlist, bins = [i for i in indarr] + [gwmat.shape[0]])
                processrowchroms.update(chrlist[np.where(rowcounts > 0)])

                # copy list to make sure the original does not get altered
                if indtype == 'exclude':
                    excluderows = excluderows.copy()

                else:
                    imputerows = imputerows.copy()

    else:
        processrowchroms = set()

    if excludecols or imputecols:
        processcolchroms = set()
        for indlist, indtype in zip([excludecols, imputecols], ['exclude', 'impute']):
            if not indlist:
                continue

            else:
                colcounts, colbins = np.histogram(indlist, bins=[i for i in indarr] + [gwmat.shape[0]])
                processcolchroms.update(chrlist[np.where(colcounts > 0)])

                # copy list to make sure the original does not get altered
                if indtype == 'exclude':
                    excludecols = excludecols.copy()

                else:
                    imputecols = imputecols.copy()

    else:
        processcolchroms = set()

    # transformed index for row and col
    rowtransform = [0]
    coltransform = [0]
    for i, chr in enumerate(chrlist):
        if chr == 'chrX' and not withX:
            continue

        else:
            if i%2 == 0:
                rowtransform.append(rowtransform[-1] + indarr[i + 1] - indarr[i] if i + 1 != len(chrlist)
                                        else rowtransform[-1] + gwmat.shape[0] - indarr[i])

                rowindex[indarr[i]: indarr[i + 1] if i + 1 != len(chrlist) else gwmat.shape[0]] = True

                if chr in processrowchroms:
                    for indlist in [imputerows, excludecols]:
                        if not indlist:
                            continue

                        else:
                            for j in range(len(indlist)):
                                if indlist[j] < indarr[i + 1] if i != len(chrlist) else gwmat.shape[0]:
                                    indlist[j] = indlist[j] - indarr[i] + rowtransform[-2]

            else:
                coltransform.append(coltransform[-1] + indarr[i + 1] - indarr[i] if i + 1 != len(chrlist)
                                        else coltransform[-1] + gwmat.shape[0] - indarr[i])

                colindex[indarr[i]: indarr[i + 1] if i + 1 != len(chrlist) else gwmat.shape[1]] = True

                if chr in processcolchroms:
                    for indlist in [imputecols, excludecols]:
                        if not indlist:
                            continue

                        else:
                            for j in range(len(indlist)):
                                if indlist[j] < indarr[i + 1] if i != len(chrlist) else gwmat.shape[0]:
                                    indlist[j] = indlist[j] - indarr[i] + coltransform[-2]

    # constructing interchromosome contact matrix
    Cij = gwmat[rowindex, :][:, colindex]

    # setting inf and nan to 0
    Cij[np.isnan(Cij) | np.isinf(Cij)] = 0

    if even:
        Cij = Cij.T

        # if even we flip excludes and imputes
        tmpexclude, tmpimpute = excluderows, imputerows
        excluderows, imputerows = excludecols, imputecols
        excludecols, imputecols = tmpexclude, tmpimpute

    # computing fractions of 0 elements in rows
    rowzerofrac = 1 - np.count_nonzero(Cij, axis = 1)/Cij.shape[1]
    colzerofrac = 1 - np.count_nonzero(Cij, axis = 0)/Cij.shape[0]

    # finding indices of rows and removing them
    rowrembins = np.where(rowzerofrac > removelim)[0]
    boolrowrembins = rowzerofrac > removelim

    if excluderows:
        rowrembins = np.concatenate([rowrembins, np.array(excluderows)])
        rowrembins.sort()
        boolrowrembins[excluderows] = True

    Cij = Cij[~boolrowrembins, :]

    # same for columns
    colrembins = np.where(colzerofrac > removelim)[0]
    boolcolrembins = colzerofrac > removelim

    if excludecols:
        colrembins = np.concatenate([colrembins, np.array(excludecols)])
        colrembins.sort()
        boolcolrembins[excludecols] = True

    Cij = Cij[:, ~boolcolrembins]

    if transform:
        # making sure logarithm is well defined
        Cij[Cij == 0] = 1

        # taking the logarithm
        Cij = np.log(Cij)

        # imputing values
        if imputerows and imputecols:
            colsubtract = (colrembins < imputecols[0]).sum()
            startcol = imputecols[0] - colsubtract
            endcol = imputecols[-1] - colsubtract - (colrembins > imputecols[0]).sum() + \
                        (colrembins > imputecols[-1]).sum() + 1

            rowsubtract = (rowrembins < imputerows[0]).sum()
            startrow = imputerows[0] - rowsubtract
            endrow = imputerows[-1] - rowsubtract - (rowrembins > imputerows[0]).sum() + \
                        (rowrembins > imputerows[-1]).sum() + 1

            for row in range(startrow, endrow):
                rv = scistats.norm(loc=Cij[row].mean(), scale=Cij[row].std())
                zeros = np.where(Cij[row] == 0)[0]
                Cij[row, startcol: endcol] = rv.rvs(size= endcol - startcol)
                Cij[row, zeros] = 0

        elif (imputerows and not imputecols) or (imputecols and not imputerows):
            raise Exception('Both imputerows and cols have to be given')

        # applying row-wise zscore calculation
        Cij = scistats.zscore(Cij, axis = 1, ddof = 1)

    return Cij, rowrembins, colrembins


def clusterMatrix(Cij, n_components, covariance_type = 'diag', n_iter = 1000):
    '''
    clusterMatrix(Cij, n_components, covariance_type = 'diag', n_iter = 1000)

    applies a GaussianHMM clustering to the processed genomewide interchromosomal
    contact matrix Cij as generated by constructClusterContactMatrix

    :param Cij:             processed genomewide interchromosomal contact matrix
                            as generated by constructClusterContactMatrix
    :param n_components:    number of clusters to find
    :param covariance_type: type of the covariance matrix to use (see hmmlearn documentation for more details)
    :param n_iter:          number of iterations allowed

    :return:                numpy.array containing numbers from 0 to n_components - 1
                            specifying the cluster to which each bin belongs and the model with
                            which it was calculated (i.e. fitted hmmlearn.hmm.GaussianHMM object)
    '''
    # initializing HMM object
    model = GaussianHMM(n_components = n_components, covariance_type = covariance_type, n_iter = n_iter, verbose = True)

    # fitting parameters
    model.fit(Cij)

    # compute the most likely state sequence using viterbi
    clusters = model.predict(Cij)

    return clusters, model


def plotClustering(Cij, clusters, ax, colors, vmin = 0, vmax = 10, title = None):
    '''
    plotClustering(Cij, clusters, ax, colors, vmin = 0, vmax = 10, title = None)

    plots Cij using the information given by clusters (i.e. row-wise clusters inferred
    by GaussianHMM) in a multicolor heatmap

    :param Cij:         processed genomewide interchromosomal contact matrix
                        as generated by constructClusterContactMatrix
    :param clusters:    numpy.array containing most probable hidden state per row in Cij
    :param ax:          ax to which the heatmap should be plotted
    :param colors:      list of colors to use for each cluster
    :param vmax:        maximum value of the array corresponding to the max color value
    :param vmin:        minimum value of the array corresponding to the min color value
    :param title:       title of the heatmap

    :return:            matplotlib.Axes object, dictionary of pcolormeshes
    '''

    meshes = {}
    for clustnum, color in zip(np.unique(clusters), colors):
        # constructing masking array
        mask = np.tile(clusters != clustnum, (Cij.shape[1], 1)).T

        # masking Cij
        maskcij = np.ma.masked_where(mask, Cij)

        # generating colormap from color
        cmap = clr.LinearSegmentedColormap.from_list('map' + str(clustnum), ['White', color], N = 256)

        # plotting heatmap
        meshes[clustnum] = ax.imshow(maskcij, cmap = cmap, label = str(clustnum), vmin = vmin, vmax = vmax, aspect = 'auto')

    if title:
       ax.set_title(title)

    return ax, meshes


def computeInformationCriteria(model, data, n, covariance_type = 'diag'):
    '''
    plotInformationCriterion(data, model, covariance_type = 'diag')

    given the original data and the model of a gaussian HMM computes the AIC and BIC

    :param model:           fitted gaussian HMM
    :param data:            data on which the gaussian HMM model was fitted
    :param n:               number of components used in model fitting (i.e. number of clusters)
    :param covariance_type: type of the covariance matrix to use (see hmmlearn documentation for more details)

    :return:                matplotlib.Axes object containing the plot, dictionary of numpy.arrays containing IC values
    '''
    # calculating number of model parameters
    # for a given model with N states the number of transition parameters
    # can be calculated by noticing that at each time t we are able to transit
    # to any other state in N. Since sum(P(S)) = 1 we know the last parameter if
    # we know the other N - 1
    transitionparams = n*(n-1)

    # each state than harbors a certain emission probability governed by
    # a multivariate gaussian distribution for which we have
    # M parameters controlling the means and in case of the diagonal
    # covariance matrix we have another M parameters (M(M + 1)/2 in case of full)
    # M is the number of variables (i.e. number of rows)
    if covariance_type == 'diag':
        emissionparams = 2*n*data.shape[0]

    else:
        emissionparams = n*data.shape[0]*(data.shape[0] + 3)/2

    BIC = -2*model.score(data) + np.log(data.shape[0])*(transitionparams + emissionparams)
    AIC = -2*model.score(data) + 2*(transitionparams + emissionparams)

    return AIC, BIC

def plotInformationCriterion(values, label, mink, maxk, title, ax):
    ax.plot(np.arange(mink, maxk + 1), values[mink - 1:], label = label, ls = '--', zorder = 1)
    ax.scatter(np.arange(mink, maxk + 1), values[mink - 1:], marker = '.', color = 'dimgrey')

    ax.set_ylabel(title)
    ax.set_xlabel('n clusters')
    ax.legend()
    ax.set_xlim(mink, maxk)
    ax.set_title(title)

    return ax


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-m', '--matrix', required = True,
                    help = '''cool, h5 or npz file holding the genomewide KR normalized contact matrix.
                              Also has to contain the indexarray and chromosome list as returned by
                              constructGenomeWideContactMatrix''')
parser.add_argument('--inputFormat', default = 'cool', choices = ['cool', 'h5', 'npz'],
                    help = 'specifies the format of the input file')
parser.add_argument('--weightName', default = 'weight',
                    help = 'name of the correction weights to apply to matrix if --inputFormat == cool')
parser.add_argument('--mink', default = 1, type = int,
                    help = 'minimum number of clusters k')
parser.add_argument('--maxk', default = 20, type = int,
                    help = 'maximumn number of clusters k')
parser.add_argument('-r', '--removelim', default = 0.3, type = float,
                    help = '''determines the fraction of entries in a row/col of the clustering matrix allowed to be 0
                              if the fraction of entries is larger the row/col is removed before clustering''')
parser.add_argument('-p', '--prefix', required = True,
                    help = '''name of the output npz file holding cluster assignments for each value of k
                              between --mink and --maxk and the correspondingly calculated BIC and AIC''')
parser.add_argument('--includeChromosome', nargs = '*', default = None,
                    help = 'chromosomes to include in the normalization as space separated list (i.e. chr1 chr2 ...)')
parser.add_argument('-e', '--exclude', default = None,
                    help = '''ranges of rows with respect to the genome wide matrix to exclude from the analysis
                              has to be passed as comma-separated list of integers delimited by a colon
                              e.g. i1:i2,i3:i4,...''')
parser.add_argument('--imputerows', default = None,
                    help = '''ranges of rows with respect to the genome wide matrix that should be imputed at --imputecols positions by row
                              normal distributed values has to be passed as integers delimited by a colon e.g. i1:i2''')
parser.add_argument('--imputecols', default = None,
                    help = '''ranges of cols with respect to the genome wide matrix that should be imputed at --imputerows positions
                              has to be passed as integers delimited by a colon e.g. i1:i2''')
parser.add_argument('-pd', '--plotdir', default = None,
                    help = 'directory to which to write the plots to')
parser.add_argument('-o', '--outdir', default = '.',
                    help = 'directory to write outputfiles to')
parser.add_argument('--noScale', default = False, action = 'store_true',
                    help = 'if set bypasses scaling of normalized matrix')
parser.add_argument('--scaleFactor', default = 100000, type = float,
                    help = 'factor used to scale the matrix')
args = parser.parse_args()

if args.plotdir == None:
    plotdir = args.outdir
else:
    plotdir = args.plotdir

imputerows, imputecols = [], []
if args.imputerows and args.imputecols:
    imputerows = list(range(*[int(i) for i in args.imputerows.split(':')]))
    imputecols = list(range(*[int(i) for i in args.imputecols.split(':')]))

elif (args.imputecols and not args.imputerows) or (args.imputerows and not args.imputecols):
    raise RuntimeError('Both, --imputerows and --imputecols are required for imputation of values')

logging.info('reading in normalized contact matrix')
if args.inputFormat == 'npz':
    npz = np.load(args.matrix)
    gwmat, indarr, chrlist = [npz[key] for key in ['cm', 'inds', 'chrlist']]

elif args.inputFormat == 'h5':
    gwmat, indarr, chrlist = loadH5(args.matrix,
                                    csr = False,
                                    includechroms = args.includeChromosome,
                                    dtype = float)

else:
    gwmat, indarr, chrlist = loadCooler(args.matrix,
                                        applyNorm = True,
                                        norm = args.weightName,
                                        includeChroms = args.includeChromosome,
                                        nans_to_zero = True)

if not args.noScale:
    gwmat *= args.scaleFactor

logging.info('constructing clustering matrices')
excluderows = []
if args.exclude:
    for r in args.exclude.split(','):
        i1, i2 = [int(i) for i in r.split(':')]
        excluderows += list(range(i1, i2))

clustmats = {'odd': None, 'even': None}
remcols = {'oddremcols': None, 'evenremcols': None}
remrows = {'oddremrows': None, 'evenremrows': None}
for key, even in zip(['even', 'odd'], [True, False]):
    clustmats[key], remrows[key + 'remrows'], remcols[key + 'remcols'] = \
        constructClusterContactMatrix(gwmat,
                                      chrlist,
                                      indarr,
                                      even = even,
                                      excluderows = excluderows.copy(),
                                      imputerows = imputerows.copy(),
                                      imputecols = imputecols.copy(),
                                      removelim = args.removelim)

    nr, nc = clustmats[key].shape
    logging.info('removed %0.2f percent of rows and %0.2f percent of cols with > %0.2f percent 0 entries for clustermatrix of %s chromosomes'
                 % ((1 - nr / (nr + len(remrows[key + 'remrows']))) * 100,
                    (1 - nc / (nc + len(remcols[key + 'remcols']))) * 100,
                    args.removelim * 100, key))


logging.info('performing clustering for clusters k between %i and %i' % (args.mink, args.maxk))
ICdict = {key: np.zeros(shape = args.maxk) for key in ['evenAIC', 'oddAIC', 'evenBIC', 'oddBIC']}
cfig, caxs = plt.subplots(1, 2)
clusterassignments = {}
models = {}
basename = ntpath.basename(args.matrix).split('.')[0]
cmap = plt.get_cmap('jet')
for clustering in ['even', 'odd']:
    for k in range(args.mink, args.maxk + 1):
        clusters, model = clusterMatrix(clustmats[clustering], k)
        ICdict[clustering + 'AIC'][k - 1], ICdict[clustering + 'BIC'][k - 1] = \
            computeInformationCriteria(model, clustmats[clustering], k)

        clusterassignments[clustering + 'k' + str(k)] = clusters
        models[clustering + 'k' + str(k)] = model

        colorlist = cmap(np.arange(k)/k)
        fig, ax = plt.subplots()
        ax, meshes = plotClustering(clustmats[clustering], clusters, ax, colorlist, vmax = 1)
        fig.set_figwidth(20)
        fig.set_figheight(20)
        fig.tight_layout()
        fig.savefig(os.path.join(plotdir, '_'.join([basename, clustering, 'k' + str(k)]) + '.pdf'))
        plt.close(fig)

    for criterion, cax in zip(['AIC', 'BIC'], caxs):
        plotInformationCriterion(ICdict[clustering + criterion],
                                 clustering,
                                 args.mink,
                                 args.maxk,
                                 criterion,
                                 cax)

cfig.set_figwidth(8)
cfig.set_figheight(4)
cfig.tight_layout()
cfig.savefig(os.path.join(args.plotdir, '_'.join([basename, 'informationcriterion.pdf'])))

logging.info('saving arrays')
np.savez(os.path.join(args.outdir, args.prefix + '.npz'),
                      **clusterassignments,
                      **ICdict,
                      **remcols,
                      **remrows)

for key, model in models.items():
    with open(os.path.join(args.outdir, args.prefix + 'hmm' + key + '.pkl'), 'wb') as file:
        pickle.dump(model, file)
