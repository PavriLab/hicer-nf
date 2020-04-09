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
import importlib
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
from scipy.sparse import triu, csr_matrix

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


class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pFileType='cool', pMatrixFile=None, pChrnameList=None,
                 pApplyCorrectionCoolerLoad=None, pBedFileHicPro=None, pCorrectionFactorTable=None,
                 pCorrectionOperator=None, pEnforceInteger=None, pAppend=None, pFileWasH5=None, pHiCInfo=None, pHic2CoolVersion=None):

        self.class_ = getattr(importlib.import_module('.' + pFileType.lower(), package='hicmatrix.lib'), pFileType.title())

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
parser = ap.ArgumentParser()
parser.add_argument('-m', '--matrix', required = True,
                    help = 'path to *.h5 matrix file generated by HiCExplorer (has to be single resolution)')
parser.add_argument('-o', '--outputFile', required = True,
                    help = 'name of the output file to write the cooler to')
args = parser.parse_args()

matrixFileHandlerInput = MatrixFileHandler(pFileType='h5', pMatrixFile=args.matrix,
                                           pCorrectionFactorTable='weight',
                                           pCorrectionOperator=None,
                                           pEnforceInteger=False,
                                           pApplyCorrectionCoolerLoad=True)

_matrix, cut_intervals, nan_bins, distance_counts, correction_factors = matrixFileHandlerInput.load()

matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool',
                                            pEnforceInteger=False,
                                            pFileWasH5=True,
                                            pHic2CoolVersion=None)

matrixFileHandlerOutput.set_matrix_variables(_matrix,
                                             cut_intervals,
                                             nan_bins,
                                             correction_factors,
                                             distance_counts)

matrixFileHandlerOutput.save(args.outputFile,
                             pSymmetric=True,
                             pApplyCorrection=False)