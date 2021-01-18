import logging
import argparse as ap
import pysam as ps


def get_header_from_samfile(file)
    f = ps.AlignmentFile(file, 'r')
    return f.header


def initialize_outputfiles(nameprefix, chromosomes, header):
    handle_dict = {}
    for chrom in chromosomes:
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

outfile_dict = initialize_outputfiles(args.outputPrefix,
                                      samheader.references,
                                      header)

for file in args.samfiles:
    samfile = ps.AlignmentFile(file, 'r')
