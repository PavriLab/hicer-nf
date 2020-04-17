import argparse as ap
import numpy as np
from scipy.sparse.csr import csr_matrix
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import tables
import sys
import logging

redmap = clr.LinearSegmentedColormap.from_list('redmap', ['White', 'Red'], N = 256)

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


def plotmatrix(mat, cmap, vmin=None, vmax=None, ax=None, xticks=None, yticks=None,
               xchroms=None, ychroms=None, title=None, remove=None, aspect = 'equal'):
    '''
    plotmatrix(mat, cmap, vmin = None, vmax = None, ax = None, xticks = None, yticks = None,
               xchroms = None, ychroms = None, title = None, remove = None)

    function to visualize a given contact matrix

    :param mat:         matrix to visualize
    :param cmap:        colormap to use
    :param vmin:        minimum value for colormap
    :param vmax:        maximum value for colormap
    :param ax:          matplotlib.Axes object to plot the heatmap to
                        if not given it is created anew
    :param xticks:      position of xticks (chromosome borders) if None xticks are disabled
    :param yticks:      position of yticks (chromosome borders) if None yticks are disabled
    :param xchroms:     names of the chromosomes on the x-axis
    :param ychroms:     names of the chromosomes on the y-axis
    :param title:       title of the plot
    :param remove:      expects a dictionary of the format {'row': numpy.array, 'col': numpy.array}
                        where the arrays hold indices of rows and columns that should be removed
                        before plotting
    :param aspect:      either equal (default) or auto, see matplotlib.pyplot.imshow for details
                        results in a higher resolution of the plot

    :return:            matplotlib.Axes and matplotlib.image objects
    '''

    if not ax:
        fig, ax = plt.subplots()

    # generating matrix with all zero diagonal for and
    # ad-hoc visualization
    zerodiag = mat.copy()
    np.fill_diagonal(zerodiag, 0)

    if not vmin:
        vmin = 0

    if not vmax:
        vmax = zerodiag.max()

    if remove:
        mat = mat.copy()

        # generating boolean indices
        remrows = np.ones(shape=mat.shape[0], dtype=bool)
        remcols = np.ones(shape=mat.shape[1], dtype=bool)

        for rowind in remove['row']:
            remrows[rowind] = False

        for colind in remove['col']:
            remcols[colind] = False

        mat = mat[remrows, :][:, remcols]

    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect = aspect)

    if title:
        ax.set_title(title)

    # deactivating default ticks
    ax.tick_params(axis='both', which='both', bottom=False, left=False,
                   labelbottom=False, labelleft=False)

    if np.any(xticks):
        ax.tick_params(axis='x', which='major', top=True, direction='in', length=2)
        ax.tick_params(axis='x', which='minor', labeltop=True)
        ax.set_xticks([i for i in xticks] + [mat.shape[1]])

    if xchroms:
        xlabelpos = []
        for i in range(len(xticks)):
            lp = (xticks[i + 1] + xticks[i]) / 2 if i + 1 != len(xticks) else (mat.shape[1] + xticks[i]) / 2
            xlabelpos.append(lp)

        ax.set_xticks(xlabelpos, minor=True)
        ax.set_xticklabels(xchroms, minor=True, fontsize=10)

    if np.any(yticks):
        ax.tick_params(axis='y', which='major', left=True, direction='in', length=2)
        ax.tick_params(axis='y', which='minor', labelleft=True)
        ax.set_yticks([i for i in yticks] + [mat.shape[0]])

    if ychroms:
        ylabelpos = []
        for i in range(len(yticks)):
            lp = (yticks[i + 1] + yticks[i]) / 2 if i + 1 != len(yticks) else (mat.shape[0] + yticks[i]) / 2
            ylabelpos.append(lp)

        ax.set_yticks(ylabelpos, minor=True)
        ax.set_yticklabels(ychroms, minor=True, fontsize=10)

    return ax, im


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-m', '--matrix', required = True,
                    help = 'matrix to plot in *.h5 format')
parser.add_argument('-c', '--chromosomes', default = None, nargs = '*',
                    help = 'space separated list of chromosomes to include in the plot')
parser.add_argument('--vMax', default = None,
                    help = 'max value of the colormap')
parser.add_argument('-o', '--out', required = True,
                    help = 'name of the plot file to be created')
args = parser.parse_args()

logging.info('reading %s' % args.matrix)
matrix, inds, chrlist = loadH5(args.matrix, includechroms = args.chromosomes, csr = False, dtype = float)

logging.info('generating plot')
fig, ax = plt.subplots()
ax, im = plotmatrix(matrix, redmap, vmax = args.vMax, ax = ax,
                    yticks = inds, xticks = inds,
                    ychroms = list(chrlist), xchroms = list(chrlist))

fig.set_figheight(14)
fig.set_figwidth(14)
fig.tight_layout()
fig.savefig(args.out)
