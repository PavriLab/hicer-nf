#!/usr/bin/env python

import argparse as ap
import pysam as ps
from io import StringIO
import pandas as pd
import os
import logging

def isBam(filename):
    return True if filename.endswith('bam') else False


def makeRecord(read1, read2):
    str1, chr1, pos1 = int(read1.is_reverse), read1.reference_name, read1.reference_start
    str2, chr2, pos2 = int(read2.is_reverse), read2.reference_name, read2.reference_start

    if chr1 < chr2:
        return ' '.join([str(str1), chr1, str(pos1), '0', str(str2), chr2, str(pos2), '1']) + '\n'

    elif chr1 > chr2:
        return ' '.join([str(str2), chr2, str(pos2), '0', str(str1), chr1, str(pos1), '1']) + '\n'

    elif chr1 == chr2:
        if pos1 < pos2:
            return ' '.join([str(str1), chr1, str(pos1), '0', str(str2), chr2, str(pos2), '1']) + '\n'

        elif pos1 > pos2:
            return ' '.join([str(str2), chr2, str(pos2), '0', str(str1), chr1, str(pos1), '1']) + '\n'

        elif pos1 == pos2:
            return ' '.join([str(str1), chr1, str(pos1), '0', str(str2), chr2, str(pos2), '1']) + '\n' \
                if str1 < str2 else ' '.join([str(str2), chr2, str(pos2), '0', str(str1), chr1, str(pos1), '1']) + '\n'


def writeChunk(buffer, filename, chunkNum):
    buffer.seek(0)
    frame = pd.read_csv(buffer,
                        sep=' ',
                        index_col = False,
                        header=None,
                        names=['str1', 'chr1', 'pos1', 'frag1', 'str2', 'chr2', 'pos2', 'frag2'])
    #print(frame.head())
    buffer = StringIO()

    frame.sort_values(by=['chr1', 'chr2', 'str1', 'str2', 'pos1'], inplace=True)
    #print(frame.head())
    frame.to_csv(filename + '.00' + str(chunkNum), sep=' ', header=False, index=False)

    return buffer


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--input', nargs = '+',
                    help = 'either complete or split HICUP file')
parser.add_argument('--bufferSize', default = 5000000, type = int,
                    help = 'number of lines to store in memory before writing to --outputFile')
parser.add_argument('--tmpDir', default = '.',
                    help = 'temporary directory to use')
parser.add_argument('-o', '--outputFilePrefix', required = True,
                    help = 'prefix of the outputfiles to write')
args = parser.parse_args()


readmode = 'rb' if isBam(args.input[0]) else 'r'
buffer = StringIO()
linesInBuffer = 0
fileNum = 0

if not os.path.exists(args.tmpDir):
    os.mkdir(args.tmpDir)

if len(args.input) == 1:
    with ps.AlignmentFile(args.input[0], readmode) as infile:
        while True:
            try:
                r1 = infile.__next__()
                r2 = infile.__next__()

            except StopIteration:
                buffer = writeChunk(buffer, os.path.join(args.tmpDir, args.outputFilePrefix), fileNum)
                linesInBuffer = 0
                fileNum += 1
                break

            #print(int(r1.is_reverse), r1.reference_name, r1.reference_start)
            #print(int(r2.is_reverse), r2.reference_name, r2.reference_start)
            #print(makeRecord(r1, r2))

            buffer.write(makeRecord(r1, r2))
            linesInBuffer += 1

            if linesInBuffer == args.bufferSize:
                buffer = writeChunk(buffer, os.path.join(args.tmpDir, args.outputFilePrefix), fileNum)
                linesInBuffer = 0
                fileNum += 1


elif len(args.input) == 2:
    infile1 = ps.AlignmentFile(args.input[0], readmode)
    infile2 = ps.AlignmentFile(args.input[1], readmode)

    while True:
        try:
            r1 = infile1.__next__()
            r2 = infile2.__next__()

        except StopIteration:
            buffer = writeChunk(buffer, os.path.join(args.tmpDir, args.outputFilePrefix), fileNum)
            linesInBuffer = 0
            fileNum += 1
            break

        #print(int(r1.is_reverse), r1.reference_name, r1.reference_start)
        #print(int(r2.is_reverse), r2.reference_name, r2.reference_start)
        #print(makeRecord(r1, r2))

        buffer.write(makeRecord(r1, r2))
        linesInBuffer += 1

        if linesInBuffer == args.bufferSize:
            buffer = writeChunk(buffer, os.path.join(args.tmpDir, args.outputFilePrefix), fileNum)
            linesInBuffer = 0
            fileNum += 1

    for file in [infile1, infile2]:
        file.close()
