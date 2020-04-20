#!/usr/bin/env python

#################################################################
#                                                               #
#   This script is a lightweight version of HiCExplorers        #
#   hicConvertFormat, works exclusively for h5 to cool          #
#   conversion and fixes a bug in the cooler conversion         #
#   causing the resulting cooler to be incompatible with        #
#   HiGlass.                                                    #
#                                                               #
#   Therefore, most of the code blow is copied from             #
#   the HiCMatrix (https://github.com/deeptools/HiCMatrix)      #
#   and HiCExplorer (https://github.com/deeptools/HiCExplorer)  #
#                                                               #
#################################################################

import os
os.environ['NUMEXPR_MAX_THREADS'] = '64'

import argparse as ap
import pandas as pd
import sys
import logging
from copy import deepcopy
from os import unlink
import math
import time
import cooler
import h5py
import tables
import numpy as np
from scipy.sparse import triu, tril, dia_matrix, csr_matrix
from collections import OrderedDict
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
from intervaltree import IntervalTree, Interval

log = logging.getLogger(__name__)


def convertNansToOnes(pArray):
    nan_elements = np.flatnonzero(np.isnan(pArray))
    if len(nan_elements) > 0:
        pArray[nan_elements] = 1.0
    return pArray


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


def toBytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    if sys.version_info[0] == 2:
        return s
    if isinstance(s, bytes):
        return s
    # if isinstance(s, np.bytes_):
    #     return np.bytes_(s)
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [toBytes(x) for x in s]
    return s


def check_chrom_str_bytes(pIteratableObj, pObj):
    # determine type
    if isinstance(pObj, list) and len(pObj) > 0:
        type_ = type(pObj[0])
    else:
        type_ = type(pObj)
    if not isinstance(type(next(iter(pIteratableObj))), type_):
        if type(next(iter(pIteratableObj))) is str:
            pObj = toString(pObj)
        elif type(next(iter(pIteratableObj))) in [bytes, np.bytes_]:
            pObj = toBytes(pObj)
    return pObj


class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pFileType='cool', pMatrixFile=None, pChrnameList=None,
                 pApplyCorrectionCoolerLoad=None, pBedFileHicPro=None, pCorrectionFactorTable=None,
                 pCorrectionOperator=None, pEnforceInteger=None, pAppend=None, pFileWasH5=None, pHiCInfo=None, pHic2CoolVersion=None):

        self.class_ = H5 if pFileType == 'h5' else Cool

        if pFileType == 'hicpro':
            self.matrixFile = self.class_(pMatrixFile=pMatrixFile, pBedFile=pBedFileHicPro)
        else:
            self.matrixFile = self.class_(pMatrixFile)
            if pFileType == 'cool':
                self.matrixFile.chrnameList = pChrnameList
                if pCorrectionFactorTable is not None:
                    self.matrixFile.correctionFactorTable = pCorrectionFactorTable
                if pCorrectionOperator is not None:
                    self.matrixFile.correctionOperator = pCorrectionOperator
                if pEnforceInteger is not None:
                    log.debug('pEnforceInteger {}'.format(pEnforceInteger))
                    self.matrixFile.enforceInteger = pEnforceInteger
                if pAppend is not None:
                    self.matrixFile.appendData = pAppend
                if pFileWasH5 is not None:
                    self.matrixFile.fileWasH5 = pFileWasH5
                log.debug('pApplyCorrectionCoolerLoad {}'.format(pApplyCorrectionCoolerLoad))
                if pApplyCorrectionCoolerLoad is not None:
                    self.matrixFile.applyCorrectionLoad = pApplyCorrectionCoolerLoad
                if pHiCInfo is not None:
                    self.hic_metadata = pHiCInfo
                log.debug('pHic2CoolVersion : {}'.format(pHic2CoolVersion))
                if pHic2CoolVersion is not None:
                    self.matrixFile.hic2cool_version = pHic2CoolVersion

    def load(self):

        return self.matrixFile.load()

    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrixFile.set_matrix_variables(pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts)

    def save(self, pName, pSymmetric, pApplyCorrection):
        self.matrixFile.save(pName, pSymmetric, pApplyCorrection)

    def load_init(self):
        pass


class hiCMatrix:
    """
    Class to handle Hi-C matrices
    contains routines to get intrachromosomal distances
    get sub matrices by chrname.
    """

    def __init__(self, pMatrixFile=None, pChrnameList=None):
        self.non_homogeneous_warning_already_printed = False
        self.bin_size = None
        self.bin_size_homogeneous = None  # track if the bins are equally spaced or not
        self.uncorrected_matrix = None

        self.matrix = None
        self.cut_intervals = None
        self.nan_bins = None
        self.correction_factors = None
        self.distance_counts = None
        # # when NaN bins are masked, this variable becomes contains the bin index
        # # needed to put the masked bins back into the matrix.
        self.orig_bin_ids = []
        self.orig_cut_intervals = []  # similar to orig_bin_ids. Used to identify the position of masked nan bins
        self.matrixFileHandler = None
        start_time = time.time()
        if pMatrixFile is not None:
            log.debug('Load self.matrixFileHandler')
            fileType = 'cool'
            if pMatrixFile.endswith('.h5'):
                fileType = 'h5'
            self.matrixFileHandler = MatrixFileHandler(pFileType=fileType, pMatrixFile=pMatrixFile, pChrnameList=pChrnameList)
            log.debug('init time: {}'.format(time.time() - start_time))
            self.matrix, self.cut_intervals, self.nan_bins, \
                self.correction_factors, self.distance_counts = self.matrixFileHandler.load()
            # if len(self.matrix.data) == 0:
            #     log.warning('No data for {}, not initialization of object. '.format(pChrnameList))
            #     self.interval_trees = None
            #     self.chrBinBoundaries = None
            #     return
            log.debug('load time: {}'.format(time.time() - start_time))
            start_time = time.time()

            log.debug('data loaded from file handler')
            if self.nan_bins is None:
                self.nan_bins = np.array([])

            self.fillLowerTriangle()
            log.debug('triangle time: {}'.format(time.time() - start_time))
            start_time = time.time()

            log.debug('fillLowerTriangle')

            self.restoreMaskedBins()
            log.debug('restoreMaskedBins: {}'.format(time.time() - start_time))
            start_time = time.time()

            log.debug('restoreMaskedBins')

            self.interval_trees, self.chrBinBoundaries = \
                self.intervalListToIntervalTree(self.cut_intervals)
            log.debug('intervalListToIntervalTree: {}'.format(time.time() - start_time))
            start_time = time.time()

            log.debug('intervalListToIntervalTree')

        elif pMatrixFile is None:
            log.debug('Only init object, no matrix given.')
        else:
            raise Exception('matrix file not given')
        log.debug('data loaded!')

    def save(self, pMatrixName, pSymmetric=True, pApplyCorrection=False, pHiCInfo=None):
        """ As an output format cooler and mcooler are supported.
        """

        if self.matrixFileHandler is None:
            fileType = 'cool'
            if pMatrixName.endswith('h5'):
                fileType = 'h5'
            self.matrixFileHandler = MatrixFileHandler(pFileType=fileType, pHiCInfo=pHiCInfo)

        self.restoreMaskedBins()
        self.matrixFileHandler.set_matrix_variables(self.matrix, self.cut_intervals, self.nan_bins,
                                                    self.correction_factors, self.distance_counts)
        if pMatrixName.endswith('cool'):
            self.matrixFileHandler.matrixFile.hic_metadata = pHiCInfo

        if pMatrixName.endswith('cool') or pMatrixName.endswith('h5'):
            self.matrixFileHandler.save(pMatrixName, pSymmetric=pSymmetric, pApplyCorrection=pApplyCorrection)

    def getInformationCoolerBinNames(self):
        log.info('The following columns are available: {}'.format(self.matrixFileHandler.matrixFile.getInformationCoolerBinNames()))

    def fillLowerTriangle(self):
        """
        checks if the matrix is complete or if only half of the matrix was saved.
        Returns a whole matrix.
        """
        # log.debug('sum of tril: {}'.format(tril(self.matrix, k=-1).sum()))
        if tril(self.matrix, k=-1).sum() == 0:
            # this case means that the lower triangle of the
            # symmetric matrix (below the main diagonal)
            # is zero. In this case, replace the lower
            # triangle using the upper triangle
            self.matrix = self.matrix + triu(self.matrix, 1).T

        # return matrix

    def setCutIntervals(self, cut_intervals):
        """
        Replace the cut_intervals of a matrix
        """

        # check that the matrix is squared
        if len(cut_intervals) != self.matrix.shape[0]:
            raise Exception("Length of cut_intervals {} does not match the "
                            "matrix size {}".format(len(cut_intervals), self.matrix.shape))

        self.cut_intervals = cut_intervals
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def setMatrix(self, matrix, cut_intervals):
        """
        Initialize a matrix with a given matrix
        and cut_intervals. Mostly useful for
        testing.
        """

        # check that the matrix is squared
        if matrix.shape[0] != matrix.shape[1]:
            raise Exception("Matrix is not squared. Shape is {}".format(matrix.shape))
        if len(cut_intervals) != matrix.shape[0]:
            raise Exception("Length of cut_intervals {} does not match the matrix size {}".format(len(cut_intervals),
                                                                                                  matrix.shape))

        self.matrix = matrix
        self.cut_intervals = cut_intervals
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def getBinSize(self):
        """
        estimates the bin size. In case the bin size
        is not equal for all bins (maybe except for the
        bin at the en of the chromosomes) a warning is issued.
        In case of uneven bins, the median is returned.
        """
        if self.bin_size is None:
            chrom, start, end, extra = zip(*self.cut_intervals)
            diff = np.array(end) - np.array(start)
            # If there is only one bin:
            if len(diff) == 1:
                self.bin_size = diff[0]
                return self.bin_size
            # If there are more bins, the diff will be compared
            # to the median of the differences between starts
            median = int(np.median(np.diff(start)))

            # check if the bin size is
            # homogeneous
            if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
                self.bin_size_homogeneous = False
                if self.non_homogeneous_warning_already_printed is False:
                    log.warning('Bin size is not homogeneous. \
                                      Median {}\n'.format(median))
                    self.non_homogeneous_warning_already_printed = True
            self.bin_size = median
        return self.bin_size

    def getMatrix(self):
        matrix = self.matrix.todense()
        if len(self.nan_bins):
            # to set NaN values the matrix type has to be
            # float. Corrected matrices are of float
            # type while uncorrected matrices are of
            # of int type
            if np.issubdtype(self.matrix, 'float') is False:
                matrix = matrix.astype(float)
            matrix[self.nan_bins, :] = np.nan
            matrix[:, self.nan_bins] = np.nan

        return matrix

    def getChrBinRange(self, chrName):
        """
        Given a chromosome name,
        This functions return the start and end bin indices in the matrix
        """

        if chrName in self.chrBinBoundaries:
            return self.chrBinBoundaries[chrName]
        else:
            raise Exception("chrName: {} not found in chrBinBoundaries"
                            "valid chromosomes are: {}"
                            .format(chrName, self.chrBinBoundaries.keys()))
            return None, None

    def getChrNames(self):
        """
        returns the names of the chromosomes
        present in the matrix
        """
        return list(self.chrBinBoundaries)

    def getBinPos(self, binIndex):
        """
        given a bin, it returns the chromosome name,
        start position and end position
        """
        if binIndex < len(self.cut_intervals):
            return self.cut_intervals[binIndex]
        else:
            raise Exception("binIndex: {} not found".format(binIndex))
            return None

    def getRegionBinRange(self, chrname, startpos, endpos):
        """
        Given a chromosome region, this function returns
        the bin indices that overlap with such region.
        """

        try:
            # chromosome_size = hic_matrix.get_chromosome_sizes()
            # chrname = check_chrom_str_bytes(self.interval_trees, chrname)
            if type(next(iter(self.interval_trees))) != type(chrname):
                if type(next(iter(self.interval_trees))) is str:
                    chrname = toString(chrname)
                elif type(next(iter(self.interval_trees))) is bytes:
                    chrname = toBytes(chrname)
                elif type(next(iter(self.interval_trees))) is np.bytes_:
                    chrname = toBytes(chrname)
            # chr_end_pos = chromosome_size[chrname]
            # self.interval_trees[chrname]
            if chrname not in self.interval_trees:
                raise Exception("chromosome: {} name not found in matrix"
                                "valid names are: {}"
                                .format(chrname, self.interval_trees.keys()))
        except KeyError as ke:
            log.exception("chromosome: {} name not found in matrix".format(chrname))
            log.exception("valid names are: ")
            log.exception(self.interval_trees.keys())
            log.exception(str(ke))

        try:
            startpos = int(startpos)
            endpos = int(endpos)
        except ValueError as ve:
            log.exception("{} or {}  are not valid "
                          "position values.".format(startpos, endpos))
            log.exception(str(ve))

        try:

            startbin = sorted(self.interval_trees[chrname][startpos:startpos + 1])[0].data
            endbin = sorted(self.interval_trees[chrname][endpos:endpos + 1])[0].data
        except IndexError:
            # log.exception("chrname: " + chrname)
            # log.exception("len intervaltree: "+len(self.interval_trees[chrname]))
            # log.exception("start and end pos:" + startpos + ":::" + endpos )
            log.exception("Index error")
            return None

        return startbin, endbin

    @staticmethod
    def getDistList(rows, cols, cut_intervals):
        """
            Given a list of rows and cols
            an array is returned containing
            the genomic distance between
            each element of the row array
            with each element of the col array.
            -1 is returned for inter-chromosomal
            interactions.
            A matching list containing the chromosome name
            is also returned
        """
        chrnamelist, startlist, endlist, extralist = zip(*cut_intervals)
        # now the distance between any two points
        # is computed and arranged such that for each
        # element of the data array, a corespondent distance is stored
        start_row = np.take(startlist, rows)
        start_col = np.take(startlist, cols)
        dist_list = start_col - start_row

        # now  all distances that are between chromosomes are removed
        # to do this I convert the array of chromosomes to
        # a array of indices. Then, when subtracting the
        # values that correspond to matrix.row and matrix.col
        # using the array of indices, any value other
        # than 0 means inter-chromosomal row,col combination.

        # chr_id_list is based on a trick using np.unique
        # to get from a list of strings
        # a list of integers
        chr_id_list = np.unique(chrnamelist, return_inverse=True)[1]

        chr_row = np.take(chr_id_list, rows)
        chr_col = np.take(chr_id_list, cols)
        chr_diff = chr_row - chr_col
        # set in dist_list array '-1' for all interchromosomal values
        dist_list[chr_diff != 0] = -1

        # make a corresponding chromosome name list
        # if filtering per chromosome is required
        chrom_list = np.take(chrnamelist, rows)
        chrom_list[chr_diff != 0] = ''

        return dist_list, chrom_list

    @staticmethod
    def fit_cut_intervals(cut_intervals):
        # check that the matrix has bins of same size
        # otherwise try to adjust the bins to
        # to match a regular binning
        if len(cut_intervals) <= 1:
            # do nothing if there is only one interval
            return cut_intervals
        chrom, start, end, extra = zip(*cut_intervals)

        median = int(np.median(np.diff(start)))
        diff = np.array(end) - np.array(start)
        # check if the bin size is homogeneous
        if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
            # set the start position of a bin to the closest multiple
            # of the median
            def snap_nearest_multiple(start_x, m):
                resi = [-1 * (start_x % m), -start_x % m]
                return start_x + resi[np.argmin(np.abs(resi))]
            start = [snap_nearest_multiple(x, median) for x in start]
            end = [snap_nearest_multiple(x, median) for x in end]
            cut_intervals = list(zip(chrom, start, end, extra))
            log.info('[getCountsByDistance] Bin size is not '
                     'homogeneous, setting \n'
                     'the bin distance to the median: {}\n'.format(median))
        return cut_intervals

    def convert_to_zscore_matrix(self, maxdepth=None, perchr=False):
        return self.convert_to_obs_exp_matrix(maxdepth=maxdepth, zscore=True, perchr=perchr)

    def convert_to_obs_exp_matrix(self, maxdepth=None, zscore=False, perchr=False):
        """
        Converts a corrected counts matrix into a
        obs / expected matrix or z-scores fast.
        The caveat is that the obs/exp or z-score are only
        computed for non-zero values, although zero values that
        are not part of the sparse matrix are considered.
        For each diagonal the mean (and std when computing z-scores) are
        calculated and then each non-zero value of the sparse matrix is
        replaced by the obs/exp or z-score.
        Parameters
        ----------
        maxdepth: maximum distance from the diagonal to consider. All contacts beyond this distance will not
                         be considered.
        zscore: if a zscore wants to be returned instead of obs/exp
        Returns
        -------
        observed / expected sparse matrix
        nans occur where the standard deviation is zero
        """

        binsize = self.getBinSize()
        max_depth_in_bins = None

        if maxdepth:
            if maxdepth < binsize:
                raise Exception("Please specify a maxDepth larger than bin size ({})".format(binsize))

            max_depth_in_bins = int(float(maxdepth * 1.5) / binsize)
            # work only with the upper matrix
            # and remove all pixels that are beyond
            # max_depth_in_bis
            # (this is done by subtracting a second sparse matrix
            # that contains only the upper matrix that wants to be removed.
            self.matrix = triu(self.matrix, k=0, format='csr') - \
                triu(self.matrix, k=max_depth_in_bins, format='csr')
        else:
            self.matrix = triu(self.matrix, k=0, format='csr')

        self.matrix.eliminate_zeros()
        depth = None
        if zscore is True:
            from scipy.sparse import diags
            m_size = self.matrix.shape[0]
            if max_depth_in_bins is not None:
                depth = max_depth_in_bins
            else:
                depth = m_size
                estimated_size_dense_matrix = m_size ** 2 * 8
                if estimated_size_dense_matrix > 100e6:
                    log.info("To compute z-scores a dense matrix is required. This will use \n"
                             "{} Mb of memory.\n To reduce memory use the maxdeph option."
                             "".format(estimated_size_dense_matrix / 1e6))

            # to compute zscore the zero values need to be accounted and the matrix
            # need to become dense. This is only practical if only up to certain distance
            # wants to be evaluated, otherwise the dense matrix is too large.
            # To make the matrix dense and keep the same computations as when
            # the matrix is sparse the following is done:
            # A sparse diagonal matrix of shape = matrix.shape is created with ones
            # (only upper triangle contains diagonals up to maxdeph)
            # This  sparse matrix is then added to self.matrix
            # then, -1 is subtracted to the self.matrix.data, thus effectively
            # adding zeros.
            diag_mat_ones = diags(np.repeat([1], m_size * depth).reshape(depth, m_size), list(range(depth)))

            self.matrix += diag_mat_ones

        from scipy.sparse import lil_matrix
        trasf_matrix = lil_matrix(self.matrix.shape)

        chr_submatrix = OrderedDict()
        cut_intervals = OrderedDict()
        chrom_sizes = OrderedDict()
        chrom_range = OrderedDict()
        if perchr:
            for chrname in self.getChrNames():
                chr_range = self.getChrBinRange(chrname)
                chr_submatrix[chrname] = self.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]].tocoo()
                cut_intervals[chrname] = [self.cut_intervals[x] for x in range(chr_range[0], chr_range[1])]
                chrom_sizes[chrname] = [chr_submatrix[chrname].shape[0]]
                chrom_range[chrname] = (chr_range[0], chr_range[1])

        else:
            chr_submatrix['all'] = self.matrix.tocoo()
            cut_intervals['all'] = self.cut_intervals
            # chrom_sizes['all'] = np.array([v[1] - v[0] for k, v in iteritems(self.chrBinBoundaries)])
            chrom_sizes['all'] = np.array([v[1] - v[0] for k, v in self.chrBinBoundaries.items()])

            chrom_range['all'] = (0, self.matrix.shape[0])

        # for chrname, submatrix in iteritems(chr_submatrix):
        for chrname, submatrix in chr_submatrix.items():

            log.info("processing chromosome {}\n".format(chrname))
            if zscore is True:
                # this step has to be done after tocoo()
                submatrix.data -= 1

            dist_list, chrom_list = self.getDistList(submatrix.row, submatrix.col,
                                                     hiCMatrix.fit_cut_intervals(cut_intervals[chrname]))

            # to get the sum of all values at a given distance I use np.bincount which
            # is quite fast. However, the input of bincount is positive integers. Moreover
            # it returns the sum for every consecutive integer, even if this is not on the list.
            # Thus, dist_list, which contains the distance in bp between any two bins is
            # converted to bin distance.

            # Because positive integers are needed we add +1 to all bin distances
            # such that the value of -1 (which means different chromosomes) can now be used

            dist_list[dist_list == -1] = -binsize
            # divide by binsize to get a list of bin distances and add +1 to remove negative values
            dist_list = (np.array(dist_list).astype(float) / binsize).astype(int) + 1

            # for each distance, return the sum of all values
            sum_counts = np.bincount(dist_list, weights=submatrix.data)
            distance_len = np.bincount(dist_list)
            # compute the average for each distance
            mat_size = submatrix.shape[0]
            mu = {}
            std = {}
            # compute mean value for each distance

            for bin_dist_plus_one, sum_value in enumerate(sum_counts):
                if maxdepth and bin_dist_plus_one == 0:  # this is for intra chromosomal counts
                    # when max depth is set, the computation
                    # of the total_intra is not accurate and is safer to
                    # output np.nan
                    mu[bin_dist_plus_one] = np.nan
                    std[bin_dist_plus_one] = np.nan
                    continue

                if bin_dist_plus_one == 0:
                    total_intra = mat_size ** 2 - sum([size ** 2 for size in chrom_sizes[chrname]])
                    diagonal_length = int(total_intra / 2)
                else:
                    # to compute the average counts per distance we take the sum_counts and divide
                    # by the number of values on the respective diagonal
                    # which is equal to the size of each chromosome - the diagonal offset (for those
                    # chromosome larger than the offset)
                    # In the following example with two chromosomes
                    # the first (main) diagonal has a size equal to the matrix (6),
                    # while the next has 1 value less for each chromosome (4) and the last one has only 2 values

                    # 0 1 2 . . .
                    # - 0 1 . . .
                    # - - 0 . . .
                    # . . . 0 1 2
                    # . . . - 0 1
                    # . . . - - 0

                    # idx - 1 because earlier the values where
                    # shifted.
                    diagonal_length = sum([size - (bin_dist_plus_one - 1) for size in chrom_sizes[chrname] if size > (bin_dist_plus_one - 1)])
                    log.debug("Type of diagonal_length {}".format(type(diagonal_length)))

                # the diagonal length should contain the number of values at a certain distance.
                # If the matrix is dense, the distance_len[bin_dist_plus_one] correctly contains the number of values
                # If the matrix is equally spaced, then, the diagonal_length as computed before is accurate.
                # But, if the matrix is both sparse and with unequal bins, then none of the above methods is
                # accurate but the the diagonal_length as computed before will be closer.
                diagonal_length = max(diagonal_length, distance_len[bin_dist_plus_one])
                log.debug("Type of diagonal_length {}".format(type(diagonal_length)))

                if diagonal_length == 0:
                    mu[bin_dist_plus_one] = np.nan
                else:
                    mu[bin_dist_plus_one] = np.float64(sum_value) / diagonal_length

                if np.isnan(sum_value):
                    log.info("nan value found for distance {}\n".format((bin_dist_plus_one - 1) * binsize))

                # if zscore is needed, compute standard deviation: std = sqrt(mean(abs(x - x.mean())**2))
                if zscore:
                    values_sqrt_diff = \
                        np.abs((submatrix.data[dist_list == bin_dist_plus_one] - mu[bin_dist_plus_one]) ** 2)
                    # the standard deviation is the sum of the differences with mu squared (value variable)
                    # plus all zeros that are not included in the sparse matrix
                    # for which the standard deviation is
                    # (0 - mu)**2 = (mu)**2
                    # The number of zeros is the diagonal length - the length of the non zero values
                    zero_values_sqrt_diff_sum = (diagonal_length - len(values_sqrt_diff)) * mu[bin_dist_plus_one] ** 2

                    _std = np.sqrt((values_sqrt_diff.sum() + zero_values_sqrt_diff_sum) / diagonal_length)
                    std[bin_dist_plus_one] = _std

            # use the expected values to compute obs/exp
            transf_ma = np.zeros(len(submatrix.data))
            for idx, value in enumerate(submatrix.data):
                if depth is not None and dist_list[idx] > depth + 1:
                    continue
                if zscore:
                    if std[dist_list[idx]] == 0:
                        transf_ma[idx] = np.nan
                    else:
                        transf_ma[idx] = (value - mu[dist_list[idx]]) / std[dist_list[idx]]
                else:
                    transf_ma[idx] = value / mu[dist_list[idx]]

            submatrix.data = transf_ma
            trasf_matrix[chrom_range[chrname][0]:chrom_range[chrname][1], chrom_range[chrname][0]:chrom_range[chrname][1]] = submatrix.tolil()

        self.matrix = trasf_matrix.tocsr()

        return self.matrix

    @staticmethod
    def dist_list_to_dict(data, dist_list):
        """
        splits data, into numeric groups defined by dist_list
        Return a dictionary containing, for
        each unique distance a dictionary
        """

        order = np.argsort(dist_list)
        dist_list = dist_list[order]
        data = data[order]

        # having the dist_list sorted, np.split
        # is used to divide the data into
        # groups that lie at the same distance, for this
        # np.diff together with np.flatnonzero is used to
        # find the indices where the distance changes.
        # the '+1' is needed because the np.diff array is
        # one element smaller than the original array, thus
        # the indices based no the np.diff array are off by 1
        # with respect to the original array
        groups = np.split(data, np.flatnonzero(np.diff(dist_list)) + 1)

        # because the dist_list is sorted
        # the order of the unique values
        # corresponds to that of the groups.
        # In other words, group[0]
        # has distance_unique[0]
        # np.sort after np.unique  in theory
        # is not needed, but just in case...
        distance_unique = np.sort(np.unique(dist_list))

        # convert to dictionary having as key
        # the distance
        distance = {}
        for index in range(len(distance_unique)):
            distance[distance_unique[index]] = groups[index]

        return distance

    def keepOnlyTheseChr(self, chromosome_list):
        """
        given a list of chromosome names,
        these are kept, while any other is removed
        from the matrix
        """
        chromosome_list = check_chrom_str_bytes(self.interval_trees, chromosome_list)

        try:
            [self.chrBinBoundaries[x] for x in chromosome_list]
        except KeyError as e:
            raise Exception("Chromosome name {} not in matrix.".format(e))

        self.restoreMaskedBins()
        size = self.matrix.shape
        # initialize a 1D array containing the columns (and rows) to
        # select. By default none are selected
        sel = np.empty(size[0], dtype=np.bool)
        sel[:] = False

        for chrName in list(self.interval_trees):
            if chrName not in chromosome_list:
                continue

            # identify start and end rows
            # of chromosomes that wants to be
            # kept
            index_start, index_end = self.getChrBinRange(chrName)
            sel[index_start:index_end] = True

        sel_id = np.flatnonzero(sel)
        mat = self.matrix[sel_id, :][:, sel_id]

        # update bin ids
        self.cut_intervals = [self.cut_intervals[x] for x in sel_id]

        # update correction factors
        if self.correction_factors is not None:
            self.correction_factors = [self.correction_factors[x] for x in sel_id]

        # keep track of nan bins
        if len(self.nan_bins):
            _temp = np.zeros(size[0])
            _temp[self.nan_bins] = 1
            _temp = _temp[sel_id]
            self.nan_bins = np.flatnonzero(_temp == 1)
        else:
            self.nan_bins = []

        self.numCols = len(sel_id)

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)
        # remove distanceCounts
        try:
            self.distance_counts = None
        except AttributeError:
            pass
        self.matrix = mat
        return self.matrix

    def diagflat(self, value=np.nan):
        """
        sets
        the matrix diagonal to np.nan
        """
        M = self.matrix.shape[0]
        diagmatrix = dia_matrix((np.repeat(value, M), 0), shape=(M, M))
        self_diag = dia_matrix(([self.matrix.diagonal()], [0]), shape=(M, M))
        # take matrix, subtract the values of the diagonal such that
        # it becomes all zeros, replace with new values by adding them
        self.matrix = self.matrix - self_diag + diagmatrix
        return self.matrix

    def filterOutInterChrCounts(self):
        """
        set all inter chromosomal counts to np.nan
        """

        ma_coo = self.matrix.tocoo()
        dist_list, _ = hiCMatrix.getDistList(ma_coo.row, ma_coo.col,
                                             self.cut_intervals)

        # set to zero all cases in which dist_list is zero
        ma_coo.data[dist_list == -1] = 0

        self.matrix = ma_coo.tocsr()
        self.matrix.eliminate_zeros()
        return self.matrix

    def setMatrixValues(self, newMatrix):
        """
        replace the current matrix values
        by the given matrix values. The
        shapes have to coincide
        """
        assert self.matrix.shape == newMatrix.shape,\
            "Given matrix has different shape. New " \
            "values need to have the same shape as previous matrix."

        self.matrix = csr_matrix(newMatrix)

    def setCorrectionFactors(self, correction_factors):
        assert len(correction_factors) == self.matrix.shape[0], \
            "length of correction factors and length of matrix are different."
        self.correction_factors = correction_factors

    def reorderChromosomes(self, new_chr_order):
        new_order = []
        new_chr_order = check_chrom_str_bytes(self.chrBinBoundaries, new_chr_order)

        for chrName in new_chr_order:
            # check that the chromosome names are valid
            if chrName not in self.chrBinBoundaries:
                raise Exception("Chromosome name '{}' not found. Please check the correct spelling "
                                "of the chromosomes and try again".format(chrName))
            orig = self.chrBinBoundaries[chrName]
            new_order.extend(list(range(orig[0], orig[1])))
        self.reorderBins(new_order)

    def reorderBins(self, new_order):
        """
        reorders the rows and colums of the
        matrix according to the new order.
        The new order can be smaller
        than the original matrix. In that
        case, the ids not in the
        new order are removed.
        """
        orig_num_rows = self.matrix.shape[0]
        self.matrix = self.matrix[new_order, :][:, new_order]
        self.cut_intervals = [self.cut_intervals[x] for x in new_order]
        # reorder the masked bins
        # keep track of nan bins
        if len(self.nan_bins):
            _temp = np.zeros(orig_num_rows)
            _temp[self.nan_bins] = 1
            _temp = _temp[new_order]
            self.nan_bins = np.flatnonzero(_temp == 1)
        else:
            self.nan_bins = []

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def maskChromosomes(self, pChromosomeList):
        mask_ids = []
        pChromosomeList = check_chrom_str_bytes(self.chrBinBoundaries, pChromosomeList)

        for chromosome in pChromosomeList:
            # check that the chromosome names are valid
            if chromosome not in self.chrBinBoundaries:
                raise Exception("Chromosome name '{}' not found. Please check the correct spelling "
                                "of the chromosomes and try again".format(chromosome))
            orig = self.chrBinBoundaries[chromosome]
            mask_ids.extend(list(range(orig[0], orig[1])))
        self.maskBins(mask_ids)

    def maskBins(self, bin_ids=None):
        """
        Mask the list of bins given. Mask means
        to remove the bins from the matrix,
        and keep the information about the intervals
        as masked
        """
        # print("self.cut_intervalsMASKBINS___START", self.cut_intervals)

        if bin_ids is None or len(bin_ids) == 0:
            return
        self.printchrtoremove(bin_ids, restore_masked_bins=False)
        try:
            # check if a masked bin already exists
            if len(self.orig_bin_ids) > 0:
                M = self.matrix.shape[0]
                previous_bin_ids = self.orig_bin_ids[M:]
                # merge new and old masked bins
                bin_ids = np.unique(np.concatenate([previous_bin_ids, self.orig_bin_ids[bin_ids]]))
                np.sort(bin_ids)
                self.restoreMaskedBins()
        except Exception:
            pass

        # join with existing nan_bins
        if self.nan_bins is not None and len(self.nan_bins) > 0:
            log.info("found existing {} nan bins that will be "
                     "included for masking ".format(len(self.nan_bins)))
            bin_ids = np.unique(np.concatenate([self.nan_bins, bin_ids]))
            self.nan_bins = []
        rows = cols = np.delete(list(range(self.matrix.shape[1])), bin_ids)

        self.matrix = self.matrix[rows, :][:, cols]

        # to keep track of removed bins
        # I add their ids to the end of the rows vector
        # to reverse the changes, I just need to do an argsort
        # to put the removed bins in place
        # log.debug("bins_ids {}".format(bin_ids))
        self.orig_bin_ids = np.concatenate([rows, bin_ids])

        new_cut_intervals = [self.cut_intervals[x] for x in rows]

        self.orig_cut_intervals = new_cut_intervals + [self.cut_intervals[x] for x in bin_ids]

        self.cut_intervals = new_cut_intervals

        self.interval_trees, self.chrBinBoundaries = self.intervalListToIntervalTree(self.cut_intervals)

        if self.correction_factors is not None:
            self.correction_factors = self.correction_factors[rows]

    def update_matrix(self, new_matrix, new_cut_intervals):
        """
        give a new matrix and list of cut intervals, the matrix,  cut intervals and
        the respective tree are updated
        :param new_matrix: now values for the sparse matrix
        :param new_cut_intervals: list of cut intervals, each entry being a tuple of the form
        (chrom, start, end, coverage)
        :return:
        """
        if len(self.orig_bin_ids) > 0:
            raise Exception("matrix contains masked bins. Restore masked bins first")

        assert len(new_cut_intervals) == new_matrix.shape[0], "matrix shape and len of cut intervals do not match"

        self.matrix = new_matrix
        self.cut_intervals = new_cut_intervals

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

        self.nan_bins = np.flatnonzero(self.matrix.sum(0).A == 0)

    def restoreMaskedBins(self):
        """
        Puts backs into the matrix the bins
        removed
        """
        if len(self.orig_bin_ids) == 0:
            return
        # the rows to add are
        # as an empty sparse matrix
        M = self.matrix.shape[0]
        N = len(self.orig_bin_ids) - M
        rows_mat = csr_matrix((N, M))
        # cols to add
        cols_mat = csr_matrix((M + N, N))

        # add the rows and cols at the end of the
        # current matrix
        self.matrix = sparse_vstack([self.matrix, rows_mat])
        self.matrix = sparse_hstack([self.matrix, cols_mat], format='csr')

        # the new matrix has the right number of cols and rows, now
        # they need to be reordered to be back in their original places
        rows = cols = np.argsort(self.orig_bin_ids)
        self.matrix = self.matrix[rows, :][:, cols]
        self.cut_intervals = [self.orig_cut_intervals[x] for x in rows]
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)
        # set as nan_bins the masked bins that were restored
        self.nan_bins = self.orig_bin_ids[M:]

        if self.correction_factors is not None:
            # add missing values as nans at end of array
            self.correction_factors = np.concatenate([self.correction_factors,
                                                      np.repeat(np.nan, N)])
            # reorder array
            self.correction_factors = self.correction_factors[rows]

        # reset orig bins ids and cut intervals
        self.orig_bin_ids = []
        self.orig_cut_intervals = []
        log.info("masked bins were restored\n")

    def reorderMatrix(self, orig, dest):
        """
        Given a matrix, a region over the diagonal is moved from
        its origin to a new destination. With this method a
        new order of the chromosomes can be produced.
        :param orig: a tuple containing the indices of the region to be moved
        :param dest: the index of the region into which to insert
                     the section moved
        """

        rows = np.delete(list(range(self.matrix.shape[1])), range(orig[0], orig[1]))

        if dest > orig[1]:
            dest = dest - (orig[1] - orig[0])

        rows = cols = np.insert(
            rows, np.repeat(dest, orig[1] - orig[0]), list(range(orig[0], orig[1])))
        self.matrix = self.matrix[rows, :][:, cols]
        self.cut_intervals = [self.cut_intervals[x] for x in rows]
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

        if self.correction_factors is not None:
            self.correction_factors = self.correction_factors[rows]
        return

    def truncTrans(self, high=0.05):
        """Truncates trans contacts to remove blowouts
        Clip high counts in trans regions (i.e. between
        chromosomes) to the max value found in the 1-high*100
        percentile
        :param:  high : float, 0<high<1, optional
            Fraction of top trans interactions to be removed
        """
        mat = self.matrix.tocoo()
        dist_list = hiCMatrix.getDistList(mat.row, mat.col, self.cut_intervals)
        if np.count_nonzero(dist_list == -1) == 0:
            return
        max_inter = np.percentile(mat.data[dist_list == -1], (100 - high))
        mat.data[(mat.data >= max_inter) & (dist_list == -1)] == max_inter

        self.setMatrixValues(mat)

    def printchrtoremove(self, to_remove, label="Number of poor regions to remove", restore_masked_bins=True):
        """
        prints out the number of bin per chromosomes
        that will be removed
        """
        cnt = {}
        try:
            self.prev_to_remove
        except Exception:
            log.debug("No self.prev_to_remove defined, defining it now.")
            self.prev_to_remove = np.array([])

        # if the same information was already printed don't
        # show it again.
        if np.array_equal(self.prev_to_remove, to_remove):
            return

        if restore_masked_bins:
            try:
                # check if a masked bin already exists
                if len(self.orig_bin_ids) > 0:
                    log.info("Masked bins already present")
                    self.restoreMaskedBins()
            except Exception:
                pass
        for idx in to_remove:
            chrom = self.cut_intervals[idx][0]
            if chrom not in cnt:
                cnt[chrom] = 0
            cnt[chrom] += 1

        log.info('{}: {} {}'.format(label, len(to_remove), cnt))
        self.prev_to_remove = to_remove

    def get_chromosome_sizes(self):
        if self.chrBinBoundaries and len(self.chrBinBoundaries) > 0:
            chrom_sizes = OrderedDict()
            # for chrom, (start_bin, end_bin) in iteritems(self.chrBinBoundaries):
            for chrom, (start_bin, end_bin) in self.chrBinBoundaries.items():

                chrom, start, end, _ = self.cut_intervals[end_bin - 1]
                chrom_sizes[chrom] = end

            return chrom_sizes

    def intervalListToIntervalTree(self, interval_list):
        """
        given an ordered list of (chromosome name, start, end)
        this is transformed to a number of interval trees,
        one for each chromosome
        """
        cut_int_tree = {}
        chrbin_boundaries = OrderedDict()
        if len(interval_list) == 0:
            log.warning("Interval list is empty")
            return cut_int_tree, chrbin_boundaries

        intval_id = 0
        chr_start_id = 0
        previous_chrom = None
        for intval in interval_list:
            chrom, start, end = intval[0:3]
            start = int(start)
            end = int(end)
            if previous_chrom != chrom:
                if previous_chrom is None:
                    previous_chrom = chrom

                chrbin_boundaries[previous_chrom] = \
                    (chr_start_id, intval_id)
                chr_start_id = intval_id
                cut_int_tree[chrom] = IntervalTree()
                previous_chrom = chrom

            cut_int_tree[chrom].add(Interval(start, end, intval_id))

            intval_id += 1
        chrbin_boundaries[chrom] = (chr_start_id, intval_id)

        return cut_int_tree, chrbin_boundaries


def check_cooler(pFileName):
    if pFileName.endswith('.cool') or cooler.io.is_cooler(pFileName) or'.mcool' in pFileName:
        return True
    return False


class MatrixFile():

    def __init__(self, pMatrixFileName=None, pBedFile=None):
        self.matrixFileName = pMatrixFileName
        self.matrix = None
        self.cut_intervals = None
        self.nan_bins = None
        self.correction_factors = None
        self.distance_counts = None
        self.bedFile = pBedFile

    def load(self):
        log.error('Not implemented')

    def save(self):
        log.error('Not implemented')

    def is_of_type(self):
        log.error('Not implemented')

    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        log.debug('Seeting matrix variables')
        self.matrix = pMatrix
        self.cut_intervals = pCutIntervals
        self.nan_bins = pNanBins
        self.correction_factors = pCorrectionFactors

        self.distance_counts = pDistanceCounts


class H5(MatrixFile, object):

    def __init__(self, pMatrixFile):
        super().__init__(pMatrixFile)

    def load(self):
        """
        Loads a matrix stored in h5 format
        :param matrix_filename:
        :return: matrix, cut_intervals, nan_bins, distance_counts, correction_factors
        """
        log.debug('Load in h5 format')

        with tables.open_file(self.matrixFileName) as f:
            parts = {}
            try:
                for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                    parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
            except Exception as e:
                log.info('No h5 file. Please check parameters concerning the file type!')
                e
            matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                                shape=parts['shape'])
            # matrix = hiCMatrix.fillLowerTriangle(matrix)
            # get intervals
            intvals = {}
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                if toString(interval_part) == toString('chr_list'):
                    chrom_list = getattr(f.root.intervals, interval_part).read()
                    intvals[interval_part] = toString(chrom_list)
                else:
                    intvals[interval_part] = getattr(f.root.intervals, interval_part).read()

            cut_intervals = list(zip(intvals['chr_list'], intvals['start_list'], intvals['end_list'], intvals['extra_list']))
            assert len(cut_intervals) == matrix.shape[0], \
                "Error loading matrix. Length of bin intervals ({}) is different than the " \
                "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

            # get nan_bins
            try:
                if hasattr(f.root, 'nan_bins'):
                    nan_bins = f.root.nan_bins.read()
                else:
                    nan_bins = np.array([])
            except Exception:
                nan_bins = np.array([])

            # get correction factors
            try:
                if hasattr(f.root, 'correction_factors'):
                    correction_factors = f.root.correction_factors.read()
                    assert len(correction_factors) == matrix.shape[0], \
                        "Error loading matrix. Length of correction factors does not" \
                        "match size of matrix"
                    correction_factors = np.array(correction_factors)
                    mask = np.isnan(correction_factors)
                    correction_factors[mask] = 0
                    mask = np.isinf(correction_factors)
                    correction_factors[mask] = 0
                else:
                    correction_factors = None
            except Exception:
                correction_factors = None

            try:
                # get correction factors
                if hasattr(f.root, 'distance_counts'):
                    distance_counts = f.root.correction_factors.read()
                else:
                    distance_counts = None
            except Exception:
                distance_counts = None
            return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    def save(self, filename, pSymmetric=True, pApplyCorrection=None):
        """
        Saves a matrix using hdf5 format
        :param filename:
        :return: None
        """
        log.debug('Save in h5 format')

        # self.restoreMaskedBins()
        if not filename.endswith(".h5"):
            filename += ".h5"

        # if the file name already exists
        # try to find a new suitable name
        if os.path.isfile(filename):
            log.warning("*WARNING* File already exists {}\n "
                        "Overwriting ...\n".format(filename))

            unlink(filename)
        if self.nan_bins is None:
            self.nan_bins = np.array([])
        elif not isinstance(self.nan_bins, np.ndarray):
            self.nan_bins = np.array(self.nan_bins)

        # save only the upper triangle of the
        if pSymmetric:
            # symmetric matrix
            matrix = triu(self.matrix, k=0, format='csr')
        else:
            matrix = self.matrix
        matrix.eliminate_zeros()

        filters = tables.Filters(complevel=5, complib='blosc')
        with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
            matrix_group = h5file.create_group("/", "matrix", )
            # save the parts of the csr matrix
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                arr = np.array(getattr(matrix, matrix_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(matrix_group, matrix_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save the matrix intervals
            intervals_group = h5file.create_group("/", "intervals", )
            chr_list, start_list, end_list, extra_list = zip(*self.cut_intervals)
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                arr = np.array(eval(interval_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(intervals_group, interval_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save nan bins
            if len(self.nan_bins):
                atom = tables.Atom.from_dtype(self.nan_bins.dtype)
                ds = h5file.create_carray(h5file.root, 'nan_bins', atom,
                                          shape=self.nan_bins.shape,
                                          filters=filters)
                ds[:] = self.nan_bins

            # save corrections factors
            if self.correction_factors is not None and len(self.correction_factors):
                self.correction_factors = np.array(self.correction_factors)
                mask = np.isnan(self.correction_factors)
                self.correction_factors[mask] = 0
                atom = tables.Atom.from_dtype(self.correction_factors.dtype)
                ds = h5file.create_carray(h5file.root, 'correction_factors', atom,
                                          shape=self.correction_factors.shape,
                                          filters=filters)
                ds[:] = np.array(self.correction_factors)

            # save distance counts
            if self.distance_counts is not None and len(self.distance_counts):
                atom = tables.Atom.from_dtype(self.distance_counts.dtype)
                ds = h5file.create_carray(h5file.root, 'distance_counts', atom,
                                          shape=self.distance_counts.shape,
                                          filters=filters)
                ds[:] = np.array(self.distance_counts)


class Cool(MatrixFile, object):

    def __init__(self, pMatrixFile=None):
        super().__init__(pMatrixFile)
        self.chrnameList = None
        self.correctionFactorTable = 'weight'
        self.correctionOperator = None
        self.enforceInteger = False
        self.appendData = False
        self.fileWasH5 = False
        self.applyCorrectionLoad = True
        # self.hic_info = {}
        self.hic_metadata = {}
        self.cool_info = None

        self.hic2cool_version = None
        self.hicmatrix_version = None
        self.scaleToOriginalRange = None
        # self.correction_factors = None

    def getInformationCoolerBinNames(self):
        return cooler.Cooler(self.matrixFileName).bins().columns.values

    def load(self):
        log.debug('Load in cool format')
        self.minValue = None
        self.maxValue = None
        if self.matrixFileName is None:
            log.info('No matrix is initialized')

        try:
            cooler_file = cooler.Cooler(self.matrixFileName)
            if 'metadata' in cooler_file.info:
                self.hic_metadata = cooler_file.info['metadata']
            else:
                self.hic_metadata = None
            self.cool_info = deepcopy(cooler_file.info)
            # log.debug("cooler_file.info {}".format(cooler_file.info))
        except Exception as e:
            log.info("Could not open cooler file. Maybe the path is wrong or the given node is not available.")
            log.info('The following file was tried to open: {}'.format(self.matrixFileName))
            log.info("The following nodes are available: {}".format(cooler.fileops.list_coolers(self.matrixFileName.split("::")[0])))
            e
        log.debug('self.chrnameList {}'.format(self.chrnameList))
        if self.chrnameList is None:
            log.debug('muh 69')

            matrixDataFrame = cooler_file.matrix(balance=False, sparse=True, as_pixels=True)
            used_dtype = np.int32
            if np.iinfo(np.int32).max < cooler_file.info['nbins']:
                used_dtype = np.int64
            count_dtype = matrixDataFrame[0]['count'].dtype
            data = np.empty(cooler_file.info['nnz'], dtype=count_dtype)
            instances = np.empty(cooler_file.info['nnz'], dtype=used_dtype)
            features = np.empty(cooler_file.info['nnz'], dtype=used_dtype)
            i = 0
            size = cooler_file.info['nbins'] // 32
            if size == 0:
                size = 1
            start_pos = 0
            while i < cooler_file.info['nbins']:
                matrixDataFrameChunk = matrixDataFrame[i:i + size]
                _data = matrixDataFrameChunk['count'].values.astype(count_dtype)
                _instances = matrixDataFrameChunk['bin1_id'].values.astype(used_dtype)
                _features = matrixDataFrameChunk['bin2_id'].values.astype(used_dtype)

                data[start_pos:start_pos + len(_data)] = _data
                instances[start_pos:start_pos + len(_instances)] = _instances
                features[start_pos:start_pos + len(_features)] = _features
                start_pos += len(_features)
                i += size
                del _data
                del _instances
                del _features

            matrix = csr_matrix((data, (instances, features)), shape=(np.int(cooler_file.info['nbins']), np.int(cooler_file.info['nbins'])), dtype=count_dtype)
            self.minValue = data.min()
            self.maxValue = data.max()

            del data
            del instances
            del features
        else:
            if len(self.chrnameList) == 1:
                try:
                    log.debug('Load data')
                    matrix = cooler_file.matrix(balance=False, sparse=True).fetch(self.chrnameList[0]).tocsr()
                    # handle the case of an empty csr matrix
                    if len(matrix.data) == 0:
                        self.minValue = 0
                        self.maxValue = 0
                    else:
                        self.minValue = matrix.data.min()
                        self.maxValue = matrix.data.max()
                except ValueError as ve:
                    log.exception("Wrong chromosome format. Please check UCSC / ensembl notation.")
                    ve
            else:
                raise Exception("Operation to load more as one region is not supported.")

        cut_intervals_data_frame = None
        correction_factors_data_frame = None

        if self.chrnameList is not None:
            if len(self.chrnameList) == 1:
                cut_intervals_data_frame = cooler_file.bins().fetch(self.chrnameList[0])

                if self.correctionFactorTable in cut_intervals_data_frame:
                    correction_factors_data_frame = cut_intervals_data_frame[self.correctionFactorTable]
            else:
                raise Exception("Operation to load more than one chr from bins is not supported.")
        else:
            if self.applyCorrectionLoad and self.correctionFactorTable in cooler_file.bins():
                correction_factors_data_frame = cooler_file.bins()[[self.correctionFactorTable]][:]

            cut_intervals_data_frame = cooler_file.bins()[['chrom', 'start', 'end']][:]

        correction_factors = None
        if correction_factors_data_frame is not None and self.applyCorrectionLoad:
            # apply correction factors to matrix
            # a_i,j = a_i,j * c_i *c_j
            matrix.eliminate_zeros()
            if len(matrix.data) > 1:

                matrix.data = matrix.data.astype(float)

                correction_factors = convertNansToOnes(np.array(correction_factors_data_frame.values).flatten())
                # apply only if there are not only 1's
                if np.sum(correction_factors) != len(correction_factors):
                    matrix.sort_indices()

                    instances, features = matrix.nonzero()
                    instances_factors = correction_factors[instances]
                    features_factors = correction_factors[features]

                    if self.correctionOperator is None:
                        if 'generated-by' in cooler_file.info:
                            log.debug('cooler_file.info[\'generated-by\'] {} {}'.format(cooler_file.info['generated-by'], type(cooler_file.info['generated-by'])))
                            generated_by = toString(cooler_file.info['generated-by'])
                            if 'hic2cool' in generated_by:

                                self.hic2cool_version = generated_by.split('-')[1]
                                if self.hic2cool_version >= '0.5':
                                    log.debug('0.5')
                                    self.correctionOperator = '/'
                                else:
                                    log.debug('0.4')

                                    self.correctionOperator = '*'
                            else:
                                self.correctionOperator = '*'

                            log.debug('hic2cool: {}'.format(self.hic2cool_version))
                            log.debug('self.correctionOperator : {}'.format(self.correctionOperator))

                            # elif 'hicmatrix' in generated_by:

                            #     self.hicmatrix_version = generated_by.split('-')[1]
                            #     if self.hicmatrix_version >= '8':
                            #         self.correctionOperator = '/'
                            #     else:
                            #         self.correctionOperator = '*'
                        else:
                            self.correctionOperator = '*'

                    instances_factors *= features_factors
                    log.debug('hic2cool: {}'.format(self.hic2cool_version))
                    log.debug('self.correctionOperator: {}'.format(self.correctionOperator))
                    if self.correctionOperator == '*':
                        matrix.data *= instances_factors
                    elif self.correctionOperator == '/':
                        matrix.data /= instances_factors

                    # if self.scaleToOriginalRange is not None:
                    min_value = matrix.data.min()
                    max_value = matrix.data.max()
                    # check if max smaller one or if not same mangnitude
                    if max_value < 1 or (np.absolute(int(math.log10(max_value)) - int(math.log10(self.maxValue))) > 1):
                        desired_range_difference = self.maxValue - self.minValue

                        min_value = matrix.data.min()
                        max_value = matrix.data.max()

                        matrix.data = (matrix.data - min_value)
                        matrix.data /= (max_value - min_value)
                        matrix.data *= desired_range_difference
                        matrix.data += self.minValue
                        self.scaleToOriginalRange = True
                        # diff_scale_factor = matrix.data.max() / max_value
                        # if self.correctionOperator == '*':
                        #     correction_factors *= diff_scale_factor
                        # if self.correctionOperator == '/':
                        #     correction_factors /= diff_scale_factor

        cut_intervals = []
        time_start = time.time()
        log.debug('Creating cut_intervals {}'.format(time_start))
        for values in cut_intervals_data_frame.values:
            cut_intervals.append(tuple([toString(values[0]), values[1], values[2], 1.0]))
        log.debug('Creating cut_intervals {} DONE'.format(time.time() - time_start))
        del cut_intervals_data_frame
        del correction_factors_data_frame
        # try to restore nan_bins.
        try:
            shape = matrix.shape[0] if matrix.shape[0] < matrix.shape[1] else matrix.shape[1]
            nan_bins = np.arange(shape)
            nan_bins = np.setdiff1d(nan_bins, matrix.indices)

        except Exception:
            nan_bins = None

        distance_counts = None

        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    def save(self, pFileName, pSymmetric=True, pApplyCorrection=True):
        log.debug('Save in cool format')

        self.matrix.eliminate_zeros()

        if self.nan_bins is not None and len(self.nan_bins) > 0 and self.fileWasH5:
            # remove nan_bins
            correction_factors = np.ones(self.matrix.shape[0])
            correction_factors[self.nan_bins] = 0
            self.matrix.sort_indices()
            _instances, _features = self.matrix.nonzero()

            instances_factors = correction_factors[_instances]
            features_factors = correction_factors[_features]

            instances_factors = np.logical_not(np.logical_or(instances_factors, features_factors))
            self.matrix.data[instances_factors] = 0
            self.matrix.eliminate_zeros()

        # set possible nans in data to 0
        mask = np.isnan(self.matrix.data)

        self.matrix.data[mask] = 0
        self.matrix.eliminate_zeros()
        # save only the upper triangle of the
        if pSymmetric:
            # symmetric matrix
            self.matrix = triu(self.matrix, format='csr')
        else:
            self.matrix = self.matrix

        self.matrix.eliminate_zeros()

        # create data frame for bins
        # self.cut_intervals is having 4 tuples, bin_data_frame should have 3.correction_factors
        # it looks like it is faster to create it with 4, and drop the last one
        # instead of handling this before.
        bins_data_frame = pd.DataFrame(self.cut_intervals, columns=['chrom', 'start', 'end', 'interactions']).drop('interactions', axis=1)
        dtype_pixel = {'bin1_id': np.int32, 'bin2_id': np.int32, 'count': np.int32}
        if self.correction_factors is not None and pApplyCorrection:
            dtype_pixel['weight'] = np.float32
            if (self.hic2cool_version is not None and self.hic2cool_version >= '0.5') or self.fileWasH5:

                log.debug('wash5 true')
                self.correction_factors = np.array(self.correction_factors).flatten()
                self.correction_factors = 1 / self.correction_factors
                mask = np.isnan(self.correction_factors)
                self.correction_factors[mask] = 0
                mask = np.isinf(self.correction_factors)
                self.correction_factors[mask] = 0
                self.correctionOperator = '*'
                log.debug('inverted correction factors')
            weight = convertNansToOnes(np.array(self.correction_factors).flatten())
            bins_data_frame = bins_data_frame.assign(weight=weight)

            log.info("Reverting correction factors on matrix...")
            instances, features = self.matrix.nonzero()
            self.correction_factors = np.array(self.correction_factors)

            # do not apply if correction factors are just 1's
            instances_factors = self.correction_factors[instances]
            features_factors = self.correction_factors[features]

            instances_factors *= features_factors

            self.matrix.data = self.matrix.data.astype(float)

            # Apply the invert operation to get the original data
            log.debug('self.correctionOperator: {}'.format(self.correctionOperator))
            log.debug('self.fileWasH5: {}'.format(self.fileWasH5))

            if self.scaleToOriginalRange:
                min_value = self.matrix.data.min()
                max_value = self.matrix.data.max()
                desired_range_difference = max_value - min_value

                self.matrix.data = (self.matrix.data - self.minValue)
                self.matrix.data /= (self.maxValue - self.minValue)
                self.matrix.data *= desired_range_difference
                self.matrix.data += min_value

            if self.correctionOperator == '*' or self.correctionOperator is None:
                self.matrix.data /= instances_factors
            elif self.correctionOperator == '/' or self.fileWasH5:
                self.matrix.data *= instances_factors

            instances_factors = None
            features_factors = None

            self.matrix.eliminate_zeros()

        log.debug('self.correction_factors {}'.format(self.correction_factors))
        log.debug('pApplyCorrection {}'.format(pApplyCorrection))

        if self.correction_factors is not None and pApplyCorrection is False:
            dtype_pixel['weight'] = np.float32
            weight = convertNansToOnes(np.array(self.correction_factors).flatten())
            bins_data_frame = bins_data_frame.assign(weight=weight)

        instances, features = self.matrix.nonzero()

        matrix_data_frame = pd.DataFrame(instances, columns=['bin1_id'], dtype=np.int32)
        del instances
        matrix_data_frame = matrix_data_frame.assign(bin2_id=features)
        del features

        if self.enforceInteger:
            dtype_pixel['count'] = np.int32
            data = np.rint(self.matrix.data)
            matrix_data_frame = matrix_data_frame.assign(count=data)
        else:
            matrix_data_frame = matrix_data_frame.assign(count=self.matrix.data)

        if not self.enforceInteger and self.matrix.dtype not in [np.int32, int]:
            #log.warning("Writing non-standard cooler matrix. Datatype of matrix['count'] is: {}".format(self.matrix.dtype))
            dtype_pixel['count'] = self.matrix.dtype
        split_factor = 1
        if len(self.matrix.data) > 1e7:
            split_factor = 1e4
            matrix_data_frame = np.array_split(matrix_data_frame, split_factor)

        if self.appendData:
            self.appendData = 'a'
        else:
            self.appendData = 'w'

        info = {}
        # these fields are created by cooler lib. Can cause errors if not deleted.
        if 'metadata' in info:
            if self.hic_metadata is None:
                self.hic_metadata = info['metadata']
            del info['metadata']
        if 'bin-size' in info:
            del info['bin-size']
        if 'bin-type' in info:
            del info['bin-type']

        info['format'] = str('HDF5::Cooler')
        info['format-url'] = str('https://github.com/mirnylab/cooler')
        info['generated-by'] = str('h5toCool.py')
        info['generated-by-cooler-lib'] = str('cooler-' + cooler.__version__)

        info['tool-url'] = str('https://github.com/deeptools/HiCMatrix')

        # info['nchroms'] = int(bins_data_frame['chrom'][:].nunique())
        # info['chromosomes'] = list(bins_data_frame['chrom'][:].unique())
        # info['nnz'] = np.string_(str(self.matrix.nnz * 2))
        # info['min-value'] = np.string_(str(matrix_data_frame['count'].min()))
        # info['max-value'] = np.string_(str(matrix_data_frame['count'].max()))
        # info['sum-elements'] = int(matrix_data_frame['count'].sum())

        if self.hic_metadata is not None and 'matrix-generated-by' in self.hic_metadata:
            info['matrix-generated-by'] = str(self.hic_metadata['matrix-generated-by'])
            del self.hic_metadata['matrix-generated-by']
        if self.hic_metadata is not None and 'matrix-generated-by-url' in self.hic_metadata:
            info['matrix-generated-by-url'] = str(self.hic_metadata['matrix-generated-by-url'])
            del self.hic_metadata['matrix-generated-by-url']
        if self.hic_metadata is not None and 'genome-assembly' in self.hic_metadata:
            info['genome-assembly'] = str(self.hic_metadata['genome-assembly'])
            del self.hic_metadata['genome-assembly']

        local_temp_dir = os.path.dirname(os.path.realpath(pFileName))
        cooler.create_cooler(cool_uri=pFileName,
                             bins=bins_data_frame,
                             pixels=matrix_data_frame,
                             mode=self.appendData,
                             dtypes=dtype_pixel,
                             ordered=True,
                             metadata=self.hic_metadata,
                             temp_dir=local_temp_dir)

        if self.appendData == 'w':
            fileName = pFileName.split('::')[0]
            with h5py.File(fileName, 'r+') as h5file:
                h5file.attrs.update(info)
                h5file.close()


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser(description='Lightweight version of hicConvertMatrix for h5 to cool or mcool conversion only')
parser.add_argument('-m', '--matrices', nargs = '+',
                    help = 'path to *.h5 matrix file generated by HiCExplorer (has to be single resolution)')
parser.add_argument('--merge', default = False, action = 'store_true',
                    help = 'instead of converting each input matrix from h5 to cool generate a multicooler file containg all matrices')
parser.add_argument('-o', '--outputFiles', nargs = '+',
                    help = 'name(s) of the output (m)cool file(s) to write (number of names has to match the number of input matrices if --merge is not set)')
args = parser.parse_args()

for i, matrix in enumerate(args.matrices):
    matrixFileHandlerInput = MatrixFileHandler(pFileType='h5', pMatrixFile=matrix,
                                               pCorrectionFactorTable='weight',
                                               pCorrectionOperator=None,
                                               pEnforceInteger=False,
                                               pApplyCorrectionCoolerLoad=True)

    _matrix, cut_intervals, nan_bins, distance_counts, correction_factors = matrixFileHandlerInput.load()

    if args.merge:
        append = False
        if i > 0:
            append = True
        hic_matrix = hiCMatrix()
        hic_matrix.setMatrix(_matrix, cut_intervals)
        bin_size = hic_matrix.getBinSize()
        matrixFileHandlerOutput = MatrixFileHandler(
            pFileType='cool', pAppend=append, pFileWasH5=True)

        matrixFileHandlerOutput.set_matrix_variables(_matrix, cut_intervals, nan_bins,
                                                     correction_factors, distance_counts)
        matrixFileHandlerOutput.save(args.outputFiles[0] + '::/resolutions/' + str(
            bin_size), pSymmetric=True, pApplyCorrection=False)


    else:
        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool',
                                                    pEnforceInteger=False,
                                                    pFileWasH5=True,
                                                    pHic2CoolVersion=None)

        matrixFileHandlerOutput.set_matrix_variables(_matrix,
                                                     cut_intervals,
                                                     nan_bins,
                                                     correction_factors,
                                                     distance_counts)

        matrixFileHandlerOutput.save(args.outputFiles[i],
                                     pSymmetric=True,
                                     pApplyCorrection=False)
