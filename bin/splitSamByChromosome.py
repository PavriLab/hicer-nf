import logging
import argparse as ap
import pysam as ps

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--samfiles', nargs = '+',
                    help = 'one or more samfiles that should be split by chromosomes')
parser.add_argument('-o', '--outputPrefix', required = True,
                    help = 'prefix of output files. can also include a path')
args = parser.parse_args()
