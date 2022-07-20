#!/usr/bin/env python

import argparse as ap
import xml.etree.ElementTree as ET
import logging

logging.basicConfig(
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)
parser = ap.ArgumentParser()
parser.add_argument(
    'xml',
    help = 'chromsize xml file as given in igenomes'
)
parser.add_argument(
    'tsv',
    help = 'chromsize tsv file'
)
args = parser.parse_args()

tree = ET.parse(args.xml)
root = tree.getroot()

with open(args.tsv, 'w') as tsv:
    for entry in root.findall('chromosome'):
        line = []
        for tag in ['contigName', 'totalBases']:
            line.append(entry.get(tag))

        tsv.write('\t'.join(line) + '\n')
