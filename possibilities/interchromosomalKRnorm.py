import argparse as ap
import numpy as np
import os, sys, gc
import cooler
from scipy.sparse import csr_matrix, triu, tril
from krbalancing import kr_balancing
import logging
import tables

def readchromlens(chromlenfile):
    '''
    readchromlens(chromlenfile)

    reads chromosome sizes of all standard chromosomes (i.e. 1 - X and Y, no random)
    from the given tab separated file which is assumed to have the following format:

    chrname     chrsize

    :param chromlenfile:    name of the tab-separated file holding chromosome sizes

    :return:                dictionary with chrname as keys and chrsize as values
    '''

    with open(chromlenfile, 'r') as chrfile:
        chromlens = {}
        for line in chrfile:
            chr, size = line.rstrip().split('\t')

            if not len(chr.split('_')) > 1:
                chromlens[chr] = int(size)

    return chromlens


def readSparseMatrixFileIntra(matfile, chromlen, resolution):
    '''
    readSparseMatrixFileIntra(matfile, chromlen, resolution)

    reads in a contact matrix for intrachromosome contacts in sparse matrix
    notation from the given file the expected file format is as follows

        i       j
    3000000 3000000 84255.0
    3000000 4000000 12458.0
    4000000 4000000 119233.0
    3000000 5000000 4240.0
    4000000 5000000 12895.0

    where each row corresponds to a given entry in the upper triangle of the contact matrix
    as follows (matrix is symmetric M_i,j = M_j,i):

    Starting from a 0-based index the first row corresponds to the entry at i = 4 and j = 4
    (i.e. forth row and column) for a contactmatrix of 1 Mb resolution. Thus, the first two
    columns describe the left boundary of the respective interacting bins and the third
    column gives the raw observed read pair counts within the respective bins.

    0M          1M          2M          3M          4M          5M
    |===========|===========|===========|===========|===========|===========
                                            --->        <---

    Therefore, the entry 3000000 4000000 12458.0 corresponds to 12458 read pairs that map
    to the bin after 3M with one end and to the bin after 4M with their other end.

    :param matfile:     name of the file containing the upper triangle of the contact matrix
                        in sparse matrix notation
    :param chromlen:    length of the corresponding chromosome
    :param resolution:  resolution of the contact matrix (i.e. size of the bins)

    :return:            numpy.array containing the symmetric contact matrix
    '''

    #initializing contact matrix
    cmat = np.zeros(shape = (np.ceil(chromlen/resolution),)*2)

    #reading file
    with open(matfile, 'r') as file:
        for line in file:
            # extracting line entries
            i, j, c = line.rstrip().split('\t')

            # converting them to original type and dividing
            # matrix indices by resolution to get correct index
            # (i = row, j = col)
            i, j, c = int(i)//resolution, int(j)//resolution, float(c)

            cmat[[i, j], [j, i]] = c

    return cmat


def readSparseMatrixFileInter(matfile, chr1len, chr2len, resolution):
    '''
    readSparseMatrixFileInter(matfile, chr1len, chr2len, resolution)

    reads in a contact matrix for interchromosome contacts in sparse matrix
    notation from the given file the expected file format is as follows

        i       j
    3000000 3000000 53.0
    4000000 3000000 89.0
    5000000 3000000 71.0
    6000000 3000000 76.0
    7000000 3000000 88.0

    where each row corresponds to a given entry in the upper triangle of the contact matrix
    as follows (matrix is symmetric M_i,j = M_j,i):

    Starting from a 0-based index the first row corresponds to the entry at i = 4 and j = 4
    (i.e. forth row and column) for a contactmatrix of 1 Mb resolution. Thus, the first two
    columns describe the left boundary of the respective interacting bins and the third
    column gives the raw observed read pair counts within the respective bins.

                0M          1M          2M          3M          4M          5M
    chrU        |===========|===========|===========|===========|===========|===========
                                                        --->

                                                          <---
    chrV        |===========|===========|===========|===========|===========|===========
                0M          1M          2M          3M          4M          5M

    Therefore, the entry 3000000 3000000 53.0 corresponds to 53 read pairs that map
    to the bin after 3M on chrU with one end and to the bin after 3M on chrV with their other end.

    :param matfile:     name of the file containing the upper triangle of the contact matrix
                        in sparse matrix notation
    :param chr1len:     length of the first chromosome (i.e. chromosome binned in column i)
    :param chr2len      length of the second chromosome (i.e. chromosome binned in column j)
    :param resolution:  resolution of the contact matrix (i.e. size of the bins)

    :return:            numpy.array containing the symmetric contact matrix
    '''

    # initializing contact matrix
    cmat = np.zeros(shape=[np.ceil(chrlen/resolution) for chrlen in [chr1len, chr2len]])

    # reading file
    with open(matfile, 'r') as file:
        for line in file:
            # extracting line entries
            i, j, c = line.rstrip().split('\t')

            # converting them to original type and dividing
            # matrix indices by resolution to get correct index
            # (i = row, j = col)
            i, j, c = int(i) // resolution, int(j) // resolution, float(c)

            cmat[i, j] = c

    return cmat


def constructGenomeWideContactMatrix(matrixdir, chromlens, resolution, includechroms = None):
    '''
    constructGenomeWideContactMatrix(matrixdir, chromlens, resolution, includechroms = None)

    parses a directory containing all intra- and interchromosomal contact matrices
    where intrachromosomal matrices are stored in a subdirectory named dir/intra
    and all interchromosomal matrices are stored in a subdirectory named dir/inter
    The function reads in all files accordingly and assembles a genomewide contact
    matrix also returning a numpy.array indicationg index boundaries of contact matrices
    sorted by chromosome name as follows (X and Y chromosomes are listed behind all other)

    chr1    chr2  ...  chr19   chrX
     i1      i2         i19     i20

    Correspondingly the returned matrix is symmetric and has the following form

            chr1    chr2    ...    chr19    chrX
    chr1     x



    chr2
      .
      .
      .
    chr19



    chrX

    Each index in the index array corresonds to the start of a new contact matrix
    i.e. i1 refers to index 0, i2 refers to the index right after the last bin of chromosome 2
    and so forth

    :param matrixdir:       directory containing interchromosomal and intrachromosomal
                            contact matrices of a given resolution
    :param chromlens:       lengths of the chromosomes as dictionary with chrX as key
                            and length in bp as value
    :param resolution:      size of the bins used to count the contacts

    :param includechroms:   chromosomes to include in matrix construction

    :return:                numpy.array containing all matrices,
                            numpy.array specifying indices of matrix boundaries
                            sorted list of chromosomes
    '''

    # reading in matrices and computing size of final matrix
    intramats, intermats = {}, {}
    matsize = 0
    for subdir in ['intra', 'inter']:
        if subdir == 'intra':
            for file in os.listdir(os.path.join(matrixdir, subdir)):
                chr = file.split('_')[0]
                if includechroms and chr not in includechroms:
                    continue

                intramats[chr] = readSparseMatrixFileIntra(os.path.join(matrixdir, subdir, file),
                                                           chromlens[chr],
                                                           resolution
                                                           )
                matsize += intramats[chr].shape[0]

        if subdir == 'inter':
            for file in os.listdir(os.path.join(matrixdir, subdir)):
                chr1, chr2 = file.split('_')[:2]
                chr2 = 'chr' + chr2

                if includechroms:
                    if chr1 not in includechroms or chr2 not in includechroms:
                        continue

                intermats['_'.join([chr1, chr2])] = readSparseMatrixFileInter(os.path.join(matrixdir, subdir, file),
                                                                              chromlens[chr1],
                                                                              chromlens[chr2],
                                                                              resolution
                                                                              )

    # initialize final matrix
    gwmat = np.zeros(shape = (matsize, )*2)

    # specifying insertion order
    intchromlist, strchromlist = [], []
    for chr in intramats.keys():
        if chr[3:].isnumeric():
            intchromlist.append(int(chr[3:]))

        else:
            strchromlist.append(chr[3:])

    chromlist = sorted(intchromlist) + sorted(strchromlist)

    # initialize index array
    # building this in advance is necessary for the
    # construction algorithm to work properly
    indarr = np.zeros(shape = len(intramats), dtype = int)
    for i, chr in enumerate(chromlist):
        if not i + 1 == len(intramats):
            indarr[i + 1] = indarr[i] + intramats['chr' + str(chr)].shape[0]

    # inserting matrices into final matrix and recording index
    for i, chr1 in enumerate(chromlist):
        intracmat = intramats['chr' + str(chr1)]

        #inserting intrachromosome contact matrix
        gwmat[indarr[i]: indarr[i] + intracmat.shape[0], indarr[i]: indarr[i] + intracmat.shape[1]] = intracmat

        # deleting matrix to free memory
        intramats[chr1] = None
        del intracmat

        for j, chr2 in zip(range(i + 1, len(intramats)), chromlist[i + 1:]):
            key = '_'.join(['chr' + str(chr1), 'chr' + str(chr2)])
            intercmat = intermats[key]

            #print(key, j, indarr[j],
            #      (indarr[i], indarr[i] + intercmat.shape[0]),
            #      (indarr[j], indarr[j] + intercmat.shape[1]),
            #      (indarr[j], indarr[j] + intercmat.shape[1]),
            #      (indarr[i], indarr[i] + intercmat.shape[0]))

            #print(intercmat)

            # inserting interchromosome contact matrix in row
            gwmat[indarr[i]: indarr[i] + intercmat.shape[0],
                  indarr[j]: indarr[j] + intercmat.shape[1]] = intercmat

            # inserting interchromosome contact matrix in column
            gwmat[indarr[j]: indarr[j] + intercmat.shape[1],
                  indarr[i]: indarr[i] + intercmat.shape[0]] = intercmat.T

            # deleting matrix to free memory
            intermats[key] = None
            del intercmat

    return gwmat, indarr, np.array(['chr' + str(i) for i in chromlist])


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

    # normvec = octave.norm_KR(gwmat[~boolremrowcols, :][:, ~boolremrowcols])
    # normgwmat = normvec.T * (gwmat[~boolremrowcols, :][:, ~boolremrowcols] * normvec)

    # expanding the resulting matrix by the removed columns
    # to make sure we don't loose bininformation
    #j = 0
    #for i in range(retgwmat.shape[0]):
        # if the current index is in the removed indices
        # we leave the row as it is
        #if i in remrowcols:
        #    continue

        # else the row was included in the normalization
        # and we replace it with normalized values
        #else:
        #    retgwmat[i, ~boolremrowcols] = normgwmat[j, :]
        #    j += 1

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


def loadCooler(cooleruri, applyNorm = False, norm = 'weight', includeChroms = None):
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

    return matrix, np.array(inds), np.array(chroms)


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()

subparser = parser.add_subparsers(title = 'modes',
                                  description = '''choose between reading matrix from separate
                                                   text files (sparse) or *.h5 file (h5)''',
                                  dest = 'mode')

sparseparser = subparser.add_parser('sparse', help = 'reads a contact matrix from separate textfiles in sparse matrix format')
h5parser = subparser.add_parser('hdf5', help = '''reads a contact matrix from file in hdf5 formats h5 (HiCExplorer output) and cooler'''
                                               '''if a cooler is used as input you have to specify the cooler uri i.e. cool::resolutions/*''')

sparseparser.add_argument('-m', '--matrixdir', required=True,
                          help='''directory containing the intra- and interchromosome contact matrices
                                  The directory has to have the following structure:

                                  matrixdir --> intra
                                            |
                                            |
                                            --> inter

                                  where intra and inter are subdirectories of matrixdir containing the raw
                                  observed contact counts for each chromosome/chromosome-pair'''
                          )
sparseparser.add_argument('-c', '--chromsizes', required=True,
                          help='''tab-separated file of the format:
                                  chrname    chrsize
                                  containing the sizes of each contained chromosome'''
                          )
sparseparser.add_argument('-r', '--resolution', required=True, type=int,
                          help='binsize with which the corresponding matrices were generated')
sparseparser.add_argument('-o', '--outfile', required=True,
                          help='file to which the normalized genomewide matrix is written in npz format')

h5parser.add_argument('-m', '--matrix', required = True,
                      help = 'contact matrix in h5 format as given by HiCExplorer')
h5parser.add_argument('--includeChromosome', nargs = '*', default = None,
                      help = 'chromosomes to include in the normalization as space separated list (i.e. chr1 chr2 ...)')
h5parser.add_argument('--inputFormat', default = 'cool', choices = ['cool', 'h5'],
                      help = 'specifies the format of the input file')
h5parser.add_argument('-o', '--outfile', required = True,
                      help = 'file to which the normalized genomewide matrix is written in npz format')

args = parser.parse_args()

logging.info('reading matrices and constructing genomewide matrix')

if args.mode == 'sparse':
    # reading chromosome sizes
    chromlens = readchromlens(args.chromsizes)
    # constructing genomewide matrix
    gwmat, indarr, chrlist = constructGenomeWideContactMatrix(args.matrixdir, chromlens, args.resolution)

else:
    if args.inputFormat == 'h5':
        gwmat, indarr, chrlist = loadH5(args.matrix, includechroms = args.includeChromosome, csr = False)

    else:
        gwmat, indarr, chrlist = loadCooler(args.matrix, includeChroms = args.includeChromosome)

logging.info('generating interchromosomal KR-normed contact matrix')
cij = getInterchromosomeContactMatrix(gwmat, indarr, chrlist)

logging.info('saving matrix to file')
np.savez(args.outfile, cij = cij, inds = indarr, chrlist = chrlist)
