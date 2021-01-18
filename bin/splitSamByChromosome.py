#!/usr/bin/env python

import logging
import argparse as ap
import pysam as ps


def get_header_from_samfile(file)
    f = ps.AlignmentFile(file, 'r')
    return f.header


def initialize_outputfiles(nameprefix, header):
    handle_dict = {}
    for chrom in header.references:
        filename = nameprefix + '_' + chrom + '.sam'
        handle_dict[chrom] = ps.AlignmentFile(filename, 'w', header = header)

    return handle_dict


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--samfiles', nargs = '+',
                    help = 'one or more samfiles that should be split by chromosomes')
parser.add_argument('-o', '--outputPrefix', required = True,
                    help = 'prefix of output files. can also include a path')
args = parser.parse_args()

samheader = get_header_from_samfile(args.samfiles[0])

outfile_handles = initialize_outputfiles(args.outputPrefix,
                                         samheader.references,
                                         header)

for file in args.samfiles:
    with ps.AlignmentFile(file, 'r') as samfile:
        while True:
            try:
                r1 = samfile.__next__()
                r2 = samfile.__next__()

            except StopIteration:
                break

            chrom1 = r1.reference_name
            chrom2 = r2.reference_name
            handle_key = chrom1 if chrom1 < chrom2 else chrom2

            for r in [r1, r2]:
                outfile_handles[handle_key].write(r)


for key, filehandle in outfile_handles.items():
    filehandle.close()
