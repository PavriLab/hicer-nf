#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2019 Tobias Neumann, Daniel Malzl
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import os
os.environ['NUMEXPR_MAX_THREADS'] = '64'

import argparse as ap
import pyBigWig as bw
import numpy as np
import scipy.stats as scistats
import scipy.linalg as scilin
from scipy.sparse import csr_matrix
import tables
import warnings
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


def readGeneTrack(bedfile, chr, chrbins, resolution):
    '''
    reads in a genetrack file in BED format and returns a numpy.array
    with 1 if a gene is in a given bin or 0 else. This function assumes
    that the bedfile is sorted by chromosome and startposition

    :param bedfile:     BED file holding gene annotation for the given genome
    :param chr:         chromosome to retrieve the genetrack for
    :param chrbins:     number of bins the chromosome is divided in
    :param resolution:  size of the bins used to construct the contact matrix

    :return:            numpy.array holding genetrack information
    '''
    # initializing return array
    genetrack = np.zeros(shape = chrbins, dtype = int)

    # reading geneTrack
    with open(bedfile, 'r') as bed:
        brake = False
        for gene in bed:
            c, start, end, name = gene.rstrip().split('\t')[:4]

            # if we are not on the right chromosome
            if not c == chr:
                # if we already parsed the desired chromosome
                if brake:
                    # we break out of the loop
                    break

                # if we not yet parsed the desired chromosome
                else:
                    # we continue
                    continue

            else:
                # computing the bins to set value 1 to
                start, end = [int(pos) for pos in [start, end]]
                i1, i2 = start//resolution, end//resolution + 1
                genetrack[i1: i2] = 1

                if not brake:
                    brake = True

    return genetrack

def correlateEigenvectorWithGeneTrack(eigenvector, genetrack):
    '''
    computed the correlation between an eigenvector of a given chromosome and
    its gene track. If the correlation is negative the sign of the eigenvector values
    is flipped.

    :param eigenvector: first eigenvector of the correlation matrix
    :param genetrack:   numpy.array holding a 1 if there is a gene in the bin and 0 otherwise

    :return:            eigenvector with flipped signs if correlation between values and genetrack
                        is negative otherwise same as given
    '''
    correlation, pval = scistats.pearsonr(eigenvector, genetrack)

    return eigenvector if correlation >= 0 else np.negative(eigenvector)

def openEigBigWigs(prefix, header, n_eigs):
    d = {}
    for i in range(1, n_eigs + 1):
        eigbw = bw.open(prefix + '.E{}.bw'.format(i), 'w')
        eigbw.addHeader(header, maxZooms = 10)
        d[i] = eigbw

    return d

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-m', '--matrix', required = True,
                    help = 'E/O transformed HIC matrix in h5 format as generated by hicTransform')
parser.add_argument('-g', '--genes', required = True,
                    help = 'BEDfile containing the gene annotation of the genome in question')
parser.add_argument('--chromLengths', required = True,
                    help = 'tab-separated file of chromosomes and their lengths')
parser.add_argument('-r', '--resolution', required = True, type = int,
                    help = 'resolution of the given HIC matrix')
parser.add_argument('-n', '--numberOfEigenvectors', default = 2, type = int,
                    help = 'number of eigenvectors to write to file. Generates a bigWig for the first n eigenvectors each')
parser.add_argument('-o', '--outputPrefix', required = True,
                    help = 'prefix for file(s) to which to write the resulting bigwig tracks to')
args = parser.parse_args()

logging.info('reading E/O matrix')
mat, indarr, chrlist = loadH5(args.matrix, dtype = float)

chrlens = []
chrlendict = {}
with open(args.chromLengths, 'r') as file:
    for line in file:
        chr, l = line.rstrip().split('\t')
        chrlens.append((chr, int(l)))
        chrlendict[chr] = int(l)

logging.info('generating bigwig(s)')
bigWigs = openEigBigWigs(args.outputPrefix, chrlens, args.numberOfEigenvectors)

for i, chr in enumerate(chrlist):
    logging.info('processing %s' % chr)
    ind1 = indarr[i]
    ind2 = indarr[i + 1] if i + 1 != len(indarr) else mat.shape[0]

    eomat = mat[ind1: ind2, ind1: ind2].toarray()
    xi, yi = np.triu_indices(eomat.shape[0], k=1)
    eomat[yi, xi] = eomat[xi, yi]

    gtrack = readGeneTrack(args.genes, chr, eomat.shape[0], args.resolution)

    # computing pearson correlation matrix
    # ignoring divide by 0 warning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pcorrmat = np.corrcoef(eomat)

    # making sure all values are well defined
    pcorrmat[np.isnan(pcorrmat)] = 0.
    pcorrmat[np.isinf(pcorrmat)] = 0.

    # computing covariance matrix
    covmat = np.cov(pcorrmat)

    # making sure all values are well defined
    covmat[np.isnan(covmat)] = 0.
    covmat[np.isinf(covmat)] = 0.

    # computing eigenvalues and eigenvectors
    lambdas, eigvs = scilin.eigh(covmat)

    # correlating first eigenvector with genetrack
    # to flip signs if correlation is negative
    for i in range(1, args.numberOfEigenvectors + 1):
        eigv = correlateEigenvectorWithGeneTrack(eigvs[:, -i], gtrack)

        if args.resolution * len(eigv) > chrlendict[chr]:
            starts = list(range(0, args.resolution * len(eigv), args.resolution))

        else:
            # include the last bin boundary as well
            starts = list(range(0, args.resolution * len(eigv) + 1, args.resolution))

        if len(starts) > len(eigv):
            starts.pop(-1)

        ends = starts.copy()
        ends.pop(0)
        ends.append(chrlendict[chr])

        bigWigs[i].addEntries([chr] * len(eigv), starts, ends = ends, values = eigv)
