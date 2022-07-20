#!/usr/bin/env python

import argparse as ap
import pysam as ps
import logging

logging.basicConfig(
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)
parser = ap.ArgumentParser()
parser.add_argument(
    'pairedSam',
    help = 'paired SAM file'
)
args = parser.parse_args()

with ps.AlignmentFile(args.pairedSam, 'r') as sam:
    for read in sam:
        print(
            read.query_name.replace('#', '$'),
            read.reference_name,
            read.reference_start + 1,
            '-' if read.is_reverse else '+',
            sep = '\t'
        )
