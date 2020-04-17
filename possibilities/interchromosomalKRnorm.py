import argparse as ap
import numpy as np
from krbalancing import kr_balancing
import tables
from scipy.sparse import csr_matrix, triu, lil_matrix
import logging

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

            if csr:
                matrix = matrix[matrixinds, :][:, matrixinds]

            else:
                matrix = matrix.toarray()
                xi, yi = np.triu_indices(matrix.shape[0], k=1)
                matrix[yi, xi] = matrix[xi, yi]

            inds = ncuts

            chr_list = filterchrs

    if not includechroms and not csr:
        x = matrix.toarray()
        xi, yi = np.triu_indices(x.shape[0], k=1)
        x[yi, xi] = x[xi, yi]
        matrix = x

    if returnintervals:
        return matrix, np.array(inds), np.array(chr_list), intervals

    else:
        return matrix, np.array(inds), np.array(chr_list)


def KRnorm(gwmat, csr = False):
    '''
    KRnorm(gwmat, csr = False)

    normalizes a given contact matrix using the algorithm decribed by Knight and Ruiz 2012
    using the original matlab code via the oct2py octave interface.
    all rows and cols having all 0 entries are removed prior to normalization and are
    filled in after normalization again to make sure we do not loose bininformation

    :param gwmat:   raw HiC contact matrix as numpy array or scipy.sparse.csr.csr_matrix
    :param csr:     if True, returned matrix is a csr matrix

    :return:        normalized HiC contact matrix
    '''

    # finding all 0 columns (equal to all 0 rows since gwmat is symmetric)
    #remrowcols = set(np.where(np.all(gwmat == 0, axis = 1))[0])
    #boolremrowcols = np.all(gwmat == 0, axis = 1)

    # performing KR normalization on 0 col/row removed gwmat
    # kr_balancing expects scipy.sparse.csr.csr_matrix
    if type(gwmat) == np.ndarray:
        logging.info('converting to sparse matrix')
        ma = csr_matrix(gwmat)

    else:
        ma = gwmat

    # initializing kr class
    logging.info('passing pointers')
    kr = kr_balancing(ma.shape[0], ma.shape[1],
                      ma.count_nonzero(), ma.indptr.astype(np.int64, copy = False),
                      ma.indices.astype(np.int64, copy = False), ma.data.astype(np.float64, copy = False))

    # computing normalization
    logging.info('start balancing')
    kr.computeKR()

    # retreiving normalized matrix
    # True is passed to rescale the normalization by
    # the sum of the original counts
    if not csr:
        normgwmat = kr.get_normalised_matrix(True).toarray()

        # since the latter only returns an upper triangle matrix
        # we fill the lower triangle with the upper transpos
        # this is allowed due to symmetry of the matrix
        xi, yi = np.triu_indices(normgwmat.shape[0], k=1)
        normgwmat[yi, xi] = normgwmat[xi, yi]

    else:
        normgwmat = kr.get_normalised_matrix(True)

    return normgwmat


def getInterchromosomeContactMatrix(gwmat, indarr, chrlist, even = False, withX = False, norm = True):
    '''
    getInterchromosomeContactMatrix(gwmat, indarr, chrlist, even = False, withX = False, norm = True)

    uses a genomewide count matrix to extract the interchromosomal contact matrix
    as defined in constructClusterContactMatrix and normalizes it using the KR algorithm

    :param rawgwmat:    observed genomewide contact matrix either normalised by KR or RAW
    :param indarr:      array containing the indices for each chromosome bins to start
    :param chrlist:     list containing ordered chromosome names
    :param even:        if True matrix is transposed to get even chromosome bins on rows
    :param withX:       if True X chromosome is also included in the output
    :param norm:        if True the input gwmat is assumed to hold the raw contact counts
                        and is normalised by interchromosomal contacts (i.e. intrachromosomal
                        contacts are removed prior to normalization) before extraction

    :return:            interchomosome normalized contact matrix
    '''

    if norm:
        # setting intrachromosomal contacts to 0
        for i in range(len(indarr)):
            if not i == len(indarr) - 1:
                gwmat[indarr[i]: indarr[i + 1], indarr[i]: indarr[i + 1]] = 0

            else:
                gwmat[indarr[i]: , indarr[i]:] = 0

        # normalizing matrix
        gwmat = KRnorm(gwmat)

    # building boolean index arrays for row and columns
    colindex = np.zeros(shape = gwmat.shape[1], dtype = bool)
    rowindex = np.zeros(shape = gwmat.shape[0], dtype = bool)

    for i, chr in enumerate(chrlist):
        if chr == 'chrX' and not withX:
            continue

        else:
            if i%2 == 0:
                rowindex[indarr[i]: indarr[i + 1] if i + 1 != len(chrlist) else gwmat.shape[0]] = True

            else:
                colindex[indarr[i]: indarr[i + 1] if i + 1 != len(chrlist) else gwmat.shape[1]] = True

    # constructing interchromosome contact matrix
    Cij = gwmat[rowindex, :][:, colindex]

    # setting inf and nan to 0
    Cij[np.isnan(Cij) | np.isinf(Cij)] = 0

    if even:
        Cij = Cij.T

    return Cij


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser(description = 'generates an unbiased KR-normalized matrix for interchromosomal contacts')
parser.add_argument('-m', '--matrix', required = True,
                    help = 'contact matrix in h5 format as given by HiCExplorer')
parser.add_argument('--includeChromosome', nargs = '*', default = None,
                    help = 'chromosomes to include in the normalization as space separated list (i.e. chr1 chr2 ...)')
parser.add_argument('-o', '--outfile', required = True,
                    help = 'file to which the normalized genomewide matrix is written in npz format')

args = parser.parse_args()

logging.info('reading matrices and constructing genomewide matrix')

gwmat, indarr, chrlist = loadH5(args.matrix, includechroms = args.includeChromosome, csr = False)

logging.info('generating interchromosomal KR-normed contact matrix')
cij = getInterchromosomeContactMatrix(gwmat, indarr, chrlist)

logging.info('saving matrix to file')
np.savez(args.outfile, cij = cij, inds = indarr, chrlist = chrlist)
