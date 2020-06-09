#!/usr/bin/env python

import argparse as ap
import pandas as pd
import pypairix as pp
import multiprocessing as mp
from functools import partial
import logging

def extract_pairix_block(pairixfilename, outputprefix, querystring):
    ifile = pp.open(pairixfilename)
    pairs = pd.DataFrame(ifile.querys2D(querystring),
                         columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2'])


    header = '## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n'
    outputsuffix = '_'.join(querystring.split('|'))
    with open(outputprefix + '_' + outputsuffix, 'w') as ofile:
        ofile.write(header)
        pairs.to_csv(ofile,
                         sep = '\t',
                         header = False,
                         index = False,
                         mode = 'a')


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('pairix',
                    help = 'pairix indexed pairs file')
parser.add_argument('outputPrefix',
                    help = 'prefix of the split files')
parser.add_argument('-p', '--nproc', default = 1, type = int,
                    help = 'number of processors to use for splitting')
args = parser.parse_args()

blocknames = pp.open(args.pairix).get_blocknames()
func = partial(extract_pairix_block,
               args.pairix,
               args.outputPrefix)

if args.nproc > 1:
    p = mp.Pool(args.nproc)
    p.imap(func,
           blocknames)
    p.close()
    p.join()

else:
    map(func, blocknames)
