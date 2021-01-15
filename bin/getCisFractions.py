#!/usr/bin/env python

import argparse as ap
import pysam as ps
import pandas as pd
import logging
import os

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--inputFile',
                    help = 'SAM file as generated by the HICUP mapper')
parser.add_argument('-r', '--reportFile',
                    help = 'file containing the deduplication report')
args = parser.parse_args()

cis_l_10 = 0
cis_g_10 = 0
with ps.AlignmentFile(args.inputFile, 'r') as inputSam:
    while True:
        try:
            read1 = inputSam.__next__()
            read2 = inputSam.__next__()
        except StopIteration:
            break

        if read1.reference_name == read2.reference_name:
            read1pos = read1.reference_start if not read1.is_reverse else read1.reference_end
            read2pos = read2.reference_start if not read2.is_reverse else read2.reference_end

            if abs(read1pos - read2pos) < 10000:
                cis_l_10 += 1

            else:
                cis_g_10 += 1

report = pd.read_csv(args.reportFile, sep = '\t')
report.loc[0, 'Cis_<10kbp_of_uniques'] = cis_l_10
report.loc[0, 'Cis_>10kbp_of_uniques'] = cis_g_10
report.to_csv(args.reportFile, sep = '\t', index = False)