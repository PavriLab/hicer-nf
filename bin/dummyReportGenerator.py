#!/usr/bin/env python

import argparse as ap
import logging

def get_read_length(fq_handle):
    name = fq_handle.readline()
    seq = fq_handle.readline()
    fq_handle.readline()
    qual = fq_handle.readline()
    return len(seq.rstrip())

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-1', '--fastq1',
                    help = 'fastq file containing reads1 to count')
parser.add_argument('-2', '--fastq2',
                    help = 'fastq file containing reads2 to count')
parser.add_argument('-o', '--outputFile',
                    help = 'name of the summaryfile to generate')
args = parser.parse_args()

readcount1 = 0
readlength1 = 0
with open(args.fastq1, 'r') as fq1:
    rl = 1
    while rl:
        readcount1 += 1
        rl = get_read_length(fq1)
        readlength1 += rl

readcount2 = 0
readlength2 = 0
with open(args.fastq2, 'r') as fq2:
    rl = 1
    while rl:
        readcount2 += 1
        rl = get_read_length(fq2)
        readlength2 += rl


readcount1 = readcount1 // 4
readcount2 = readcount2 // 4
header = '\t'.join(['File', 'Total_Reads_Processed', 'Truncated',
                    '%Truncated', 'Not_truncated', '%Not_truncated',
                    'Average_length_truncated_sequence'])

with open(args.outputFile, 'w') as o:
    o.write(header + '\n')
    for fq, rc, rl in [(args.fastq1, readcount1, readlength1),
                       (args.fastq2, readcount2, readlength2)]:
        o.write(fq)
        for data in [rc, 0, 0, rc, 100, rl/rc]:
            o.write('\t' + str(data))

        o.write('\n')
