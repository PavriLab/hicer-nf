#!/usr/bin/env python

import argparse as ap
import pandas as pd
import logging
import datetime
from collections import defaultdict
from os import path

def initialize_frame(processedFilePrefix, reportType):
    if reportType == 'truncater':
        header = [
            'Total_Reads_Processed',
            'Truncated',
            '%Truncated',
            'Not_truncated',
            '%Not_truncated',
            'Average_length_truncated_sequence'
        ]

        df = pd.DataFrame(
            [[0, 0, .0, 0, .0, .0] for i in [1, 2]],
            columns = header,
            index = pd.Index(
                [processedFilePrefix + f'_trimmed_val_{i}.fq.gz' for i in [1, 2]],
                name = 'File'
            )
        )

    elif reportType == 'mapper':
        header = [
            'Total_reads_processed',
            'Reads_too_short_to_map',
            '%Reads_too_short_to_map',
            'Unique_alignments',
            '%Unique_alignments',
            'Multiple_alignments',
            '%Multiple_alignments',
            'Failed_to_align',
            '%failed_to_align',
            'Paired',
            '%Paired'
        ]

        df = pd.DataFrame(
            [[0, 0, .0, 0, .0, 0, .0, 0, .0, 0, .0] for i in [1, 2]],
            columns = header,
            index = pd.Index(
                [processedFilePrefix + f'_{i}.trunc.fastq' for i in [1, 2]],
                name = 'File'
            )
        )

    elif reportType == 'filter':
        header = [
            'Total_pairs',
            'Valid_pairs',
            'Cis_<10kbp',
            'Cis_>10kbp',
            'Trans',
            'Invalid_pairs',
            'Same_circularised',
            'Same_dangling_ends',
            'Same_internal',
            'Re-ligation',
            'Contiguous_sequence',
            'Wrong_size'
        ]

        df = pd.DataFrame(
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in [1, 2]],
            columns = header,
            index = pd.Index(
                [processedFilePrefix + '_1_2.pair.sam', 'dummy'],
                name = 'File'
            )
        )

    elif reportType == 'deduplicator':
        header = [
            'Read_pairs_processed',
            'Unique_di-tags',
            'Cis_<10kbp_of_uniques',
            'Cis_>10kbp_of_uniques',
            'Trans_of_uniques'
        ]

        df = pd.DataFrame(
            [[0, 0, 0, 0, 0] for i in [1, 2]],
            columns = header,
            index = pd.Index(
                [processedFilePrefix + '_1_2.filt.sam', 'dummy'],
                name = 'File'
            )
        )

    elif reportType == 'ditag_length':
        df = defaultdict(int)

    else:
        raise Exception(
            'Please specify a valid report type. Valid types are: truncater, mapper, filter, deduplicator, ditag_length'
        )

    return df


def sum_frames(df1, df2):
    for idx, row in df1.iterrows():
        if idx == 'dummy':
            continue

        for column in df1.columns:
            if column == 'File':
                continue

            else:
                df1.loc[idx, column] += df2.loc[idx, column]

    return df1


def sum_dicts(d1, d2):
    for k, v in d2.items():
        d1[k] += v

    return d1


def series_to_tsv(outputFileName, series):
    columns = ['File'] + list(series.keys())
    values = list(series.to_frame().transpose().index) + [str(int(i)) for i in series.values]
    with open(outputFileName, 'w') as ofile:
        ofile.write('\t'.join(columns) + '\n')
        ofile.write('\t'.join(values) + '\n')


def input_ditag_length_data_to_string(input_ditag_length_data):
    return ','.join([str(i) for i in input_ditag_length_data])



def generate_html_report(d, htmlstring, outputFile):
    for k, v in d.items():
        if k == 'INPUT_DITAG_LENGTH_DATA':
            v = input_ditag_length_data_to_string(v)

        htmlstring = htmlstring.replace(k, str(v))

    with open(outputFile, 'w') as htmlReport:
        htmlReport.write(htmlstring)


def get_html_template(templateFile):
    with open(templateFile, 'r') as template:
        htmlstring = template.read()

    return htmlstring


def write_summary_table(d, inputFileName, summaryFileName):
    header = [
        'File',
        'Total_Reads_1',
        'Total_Reads_2',
        'Not_Truncated_Reads_1',
        'Not_Truncated_Reads_2',
        'Truncated_Read_1',
        'Truncated_Read_2',
        'Average_Length_Truncated_1',
        'Average_Length_Truncated_2',
        'Too_Short_To_Map_Read_1',
        'Too_Short_To_Map_Read_2',
        'Unique_Alignments_Read_1',
        'Unique_Alignments_Read_2',
        'Multiple_Alignments_Read_1',
        'Multiple_Alignments_Read_2',
        'Failed_To_Align_Read_1',
        'Failed_To_Align_Read_2',
        'Paired_Read_1',
        'Paired_Read_2',
        'Valid_Pairs',
        'Invalid_Pairs',
        'Same_Circularised',
        'Same_Dangling_Ends',
        'Same_Fragment_Internal',
        'Re_Ligation',
        'Contiguous_Sequence',
        'Wrong_Size',
        'Deduplication_Read_Pairs_Uniques',
        'Deduplication_Cis_Close_Uniques',
        'Deduplication_Cis_Far_Uniques',
        'Deduplication_Trans_Uniques',
        'Percentage_Mapped',
        'Percentage_Valid',
        'Percentage_Uniques',
        'Percentage_Unique_Trans',
        'Percentage_Ditags_Passed_Through_HiCUP'
    ]

    keys = [
        'INPUT_TOTAL_READS_1',
        'INPUT_TOTAL_READS_2',
        'INPUT_NOT_TRUNCATED_READ1',
        'INPUT_NOT_TRUNCATED_READ2',
        'INPUT_TRUNCATED_READ1',
        'INPUT_TRUNCATED_READ2',
        'INPUT_AVERAGE_LENGTH_TRUNCATED_READ1',
        'INPUT_AVERAGE_LENGTH_TRUNCATED_READ2',
        'INPUT_TOO_SHORT_TO_MAP_READ_1',
        'INPUT_TOO_SHORT_TO_MAP_READ_2',
        'INPUT_UNIQUE_ALIGNMENTS_READ1',
        'INPUT_UNIQUE_ALIGNMENTS_READ2',
        'INPUT_MULTIPLE_ALIGNMENTS_READ1',
        'INPUT_MULTIPLE_ALIGNMENTS_READ2',
        'INPUT_FAILED_TO_ALIGN_READ1',
        'INPUT_FAILED_TO_ALIGN_READ2',
        'INPUT_PAIRED_READ1',
        'INPUT_PAIRED_READ2',
        'INPUT_VALID_PAIRS',
        'INPUT_INVALID_PAIRS',
        'INPUT_SAME_CIRCULARISED',
        'INPUT_SAME_DANGLING_ENDS',
        'INPUT_SAME_FRAGMENT_INTERNAL',
        'INPUT_RE_LIGATION',
        'INPUT_CONTIGUOUS_SEQUENCE',
        'INPUT_WRONG_SIZE',
        'INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES',
        'INPUT_DEDUPLICATION_CIS_CLOSE_UNIQUES',
        'INPUT_DEDUPLICATION_CIS_FAR_UNIQUES',
        'INPUT_DEDUPLICATION_TRANS_UNIQUES'
    ]

    values = [inputFileName] + \
             [str(d[k]) for k in keys] + \
             [
                str(round(d['INPUT_TOTAL_PAIRS'] / d['INPUT_TOTAL_READS_1'] * 100, 2)),
                str(round(d['INPUT_VALID_PAIRS'] / d['INPUT_TOTAL_PAIRS'] * 100, 2)),
                str(round(d['INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES'] / d['INPUT_DEDUPLICATION_READ_PAIRS_ALL'] * 100, 2)),
                str(round(d['INPUT_DEDUPLICATION_TRANS_UNIQUES'] / d['INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES'] * 100, 2)),
                str(round(d['INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES'] / d['INPUT_TOTAL_READS_1'] * 100, 2))
             ]

    with open(summaryFileName, 'w') as summaryFile:
        for line in [header, values]:
            summaryFile.write('\t'.join(line) + '\n')


logging.basicConfig(
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)
parser = ap.ArgumentParser()
parser.add_argument(
    'reports',
    nargs = '+',
    help = 'hicup report files to merge'
)
parser.add_argument(
    'template',
    help = 'hicup report html template'
)
parser.add_argument(
    'prefix',
    help = 'prefix of the merged reports'
)
parser.add_argument(
    '--writeSuppTables',
    default = False,
    action = 'store_true',
    help = 'if set also writes merged stats for each hicup step'
)
parser.add_argument(
    '-o', '--outputDir',
    default = '.',
    help = 'optional outputdirectory for the generated files'
)
args = parser.parse_args()

reportTypes = ['truncater', 'mapper', 'filter', 'deduplicator', 'ditag_length']
fileNums = {reportType: 0 for reportType in reportTypes}
statFrames = {reportType: initialize_frame(args.prefix, reportType) for reportType in reportTypes}

for reportFileName in args.reports:
    if 'ditag_size_distribution' in reportFileName:
        with open(reportFileName, 'r') as reportFile:
            d = {}
            for line in reportFile:
                k, v = [int(i) for i in line.rstrip().split('\t')]
                d[k] = v

        ditag_size_dict = sum_dicts(statFrames['ditag_length'], d)

    else:
        df = pd.read_csv(
            reportFileName,
            sep = '\t'
        )

        # making sure NaNs do not mess up the report
        # NaNs might be present in samples where there are very few reads
        # so usually not a problem
        df.fillna(0, inplace = True)
        
        df.set_index(
            'File',
            inplace = True
        )

        if 'truncater' in reportFileName:
            idx = [args.prefix + f'_trimmed_val_{i[-4: -3]}.fq.gz' for i in df.index]
            df.index = idx
            statFrames['truncater'] = sum_frames(statFrames['truncater'], df)
            fileNums['truncater'] += 1

        elif 'mapper' in reportFileName:
            idx = [args.prefix + f'_{i[-13: -12]}.trunc.fastq' for i in df.index]
            df.index = idx
            statFrames['mapper'] = sum_frames(statFrames['mapper'], df)
            fileNums['mapper'] += 1

        elif 'filter' in reportFileName:
            idx = [args.prefix + '_1_2.pair.sam']
            df.index = idx
            statFrames['filter'] = sum_frames(statFrames['filter'], df)
            fileNums['filter'] += 1

        elif 'deduplicator' in reportFileName:
            if not df.empty:
                idx = [args.prefix + '_1_2.filt.sam']
                df.index = idx
                statFrames['deduplicator'] = sum_frames(statFrames['deduplicator'], df)
                fileNums['deduplicator'] += 1


# recompute percentages
trunc_df = statFrames['truncater']
for idx, row in trunc_df.iterrows():
    for column in ['Truncated', 'Not_truncated']:
        trunc_df.loc[idx, f'%{column}'] = round(row[column] / row['Total_Reads_Processed'] * 100, 2)
    trunc_df.loc[idx, 'Average_length_truncated_sequence'] = round(row['Average_length_truncated_sequence'] / fileNums['truncater'], 2)

statFrames['truncater'] = trunc_df

map_df = statFrames['mapper']
for idx, row in map_df.iterrows():
    if idx == 'dummy':
        continue

    for column in ['Reads_too_short_to_map', 'Unique_alignments', 'Multiple_alignments', 'Failed_to_align', 'Paired']:
        if column == 'Failed_to_align':
            map_df.loc[idx, f'%{column.lower()}'] = round(row[column] / row['Total_reads_processed'] * 100, 2)

        else:
            map_df.loc[idx, f'%{column}'] = round(row[column] / row['Total_reads_processed'] * 100, 2)

statFrames['mapper'] = map_df

d = dict(
    INPUT_FILENAME = args.prefix,
    INPUT_TOTAL_READS_1 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_1.fq.gz', 'Total_Reads_Processed'],
    INPUT_TOTAL_READS_2 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_2.fq.gz', 'Total_Reads_Processed'],
    INPUT_NOT_TRUNCATED_READ1 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_1.fq.gz', 'Not_truncated'],
    INPUT_TRUNCATED_READ1 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_1.fq.gz', 'Truncated'],
    INPUT_TOO_SHORT_TO_MAP_READ_1 = statFrames['mapper'].loc[args.prefix + '_1.trunc.fastq', 'Reads_too_short_to_map'],
    INPUT_UNIQUE_ALIGNMENTS_READ1 = statFrames['mapper'].loc[args.prefix + '_1.trunc.fastq', 'Unique_alignments'],
    INPUT_MULTIPLE_ALIGNMENTS_READ1 = statFrames['mapper'].loc[args.prefix + '_1.trunc.fastq', 'Multiple_alignments'],
    INPUT_FAILED_TO_ALIGN_READ1 = statFrames['mapper'].loc[args.prefix + '_1.trunc.fastq', 'Failed_to_align'],
    INPUT_PAIRED_READ1 = statFrames['mapper'].loc[args.prefix + '_1.trunc.fastq', 'Paired'],
    INPUT_NOT_TRUNCATED_READ2 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_2.fq.gz', 'Not_truncated'],
    INPUT_TRUNCATED_READ2 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_2.fq.gz', 'Truncated'],
    INPUT_TOO_SHORT_TO_MAP_READ_2 = statFrames['mapper'].loc[args.prefix + '_2.trunc.fastq', 'Reads_too_short_to_map'],
    INPUT_UNIQUE_ALIGNMENTS_READ2 = statFrames['mapper'].loc[args.prefix + '_2.trunc.fastq', 'Unique_alignments'],
    INPUT_MULTIPLE_ALIGNMENTS_READ2 = statFrames['mapper'].loc[args.prefix + '_2.trunc.fastq', 'Multiple_alignments'],
    INPUT_FAILED_TO_ALIGN_READ2 = statFrames['mapper'].loc[args.prefix + '_2.trunc.fastq', 'Failed_to_align'],
    INPUT_PAIRED_READ2 = statFrames['mapper'].loc[args.prefix + '_2.trunc.fastq', 'Paired'],
    INPUT_TOTAL_PAIRS = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Total_pairs'],
    INPUT_VALID_PAIRS = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Valid_pairs'],
    INPUT_INVALID_PAIRS = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Invalid_pairs'],
    INPUT_SAME_CIRCULARISED = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Same_circularised'],
    INPUT_SAME_DANGLING_ENDS = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Same_dangling_ends'],
    INPUT_SAME_FRAGMENT_INTERNAL = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Same_internal'],
    INPUT_RE_LIGATION = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Re-ligation'],
    INPUT_CONTIGUOUS_SEQUENCE = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Contiguous_sequence'],
    INPUT_WRONG_SIZE = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Wrong_size'],
    INPUT_WIDTH_DITAG_LENGTH_SHORTEST = 0,
    INPUT_DITAG_LENGTH_SHORTEST = 0,
    INPUT_WIDTH_DITAG_LENGTH_LONGEST = 0,
    INPUT_DITAG_LENGTH_LONGEST = 0,
    INPUT_DITAG_LENGTH_DATA = sorted([[k,v] for k,v in statFrames['ditag_length'].items()], key = lambda x: x[0]),
    INPUT_DEDUPLICATION_READ_PAIRS_ALL = statFrames['deduplicator'].loc[args.prefix + '_1_2.filt.sam', 'Read_pairs_processed'],
    INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES = statFrames['deduplicator'].loc[args.prefix + '_1_2.filt.sam', 'Unique_di-tags'],
    INPUT_DEDUPLICATION_CIS_CLOSE_ALL = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Cis_<10kbp'],
    INPUT_DEDUPLICATION_CIS_CLOSE_UNIQUES = statFrames['deduplicator'].loc[args.prefix + '_1_2.filt.sam', 'Cis_<10kbp_of_uniques'],
    INPUT_DEDUPLICATION_CIS_FAR_ALL = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Cis_>10kbp'],
    INPUT_DEDUPLICATION_CIS_FAR_UNIQUES = statFrames['deduplicator'].loc[args.prefix + '_1_2.filt.sam', 'Cis_>10kbp_of_uniques'],
    INPUT_DEDUPLICATION_TRANS_ALL = statFrames['filter'].loc[args.prefix + '_1_2.pair.sam', 'Trans'],
    INPUT_DEDUPLICATION_TRANS_UNIQUES = statFrames['deduplicator'].loc[args.prefix + '_1_2.filt.sam', 'Trans_of_uniques'],
    INPUT_AVERAGE_LENGTH_TRUNCATED_READ1 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_1.fq.gz', 'Average_length_truncated_sequence'],
    INPUT_AVERAGE_LENGTH_TRUNCATED_READ2 = statFrames['truncater'].loc[args.prefix + '_trimmed_val_2.fq.gz', 'Average_length_truncated_sequence'],
    DATESTAMP = str(datetime.datetime.now())
)

d['INPUT_DEDUPLICATION_PERCENTAGE_UNIQUES'] = round(
    d['INPUT_DEDUPLICATION_READ_PAIRS_UNIQUES']/d['INPUT_DEDUPLICATION_READ_PAIRS_ALL'] * 100, 2
)

generate_html_report(
    d,
    get_html_template(args.template),
    path.join(args.outputDir, args.prefix + '.HiCUP_summary_report.html')
)

write_summary_table(
    d,
    args.prefix + '.hicup.sam',
    path.join(args.outputDir, f'HiCUP_summary_report_{args.prefix}.txt')
)


if args.writeSuppTables:
    for reportType, frame in statFrames.items():
        summaryFileName = path.join(
            args.outputDir,
            '_'.join([args.prefix, reportType, 'summary.txt'])
        )

        if reportType == 'filter':
            series = frame.loc[args.prefix + '_1_2.pair.sam', :]
            series_to_tsv(summaryFileName, series)

        elif reportType == 'ditag_length':
            summaryFileName = path.join(
                args.outputDir,
                '.'.join([args.prefix, 'ditag_size_distribution'])
            )
            with open(summaryFileName, 'w') as ofile:
                for k, v in frame.items():
                    ofile.write('\t'.join([str(k), str(v)]) + '\n')

        elif reportType == 'deduplicator':
            series = frame.loc[args.prefix + '_1_2.filt.sam', :]
            series_to_tsv(summaryFileName, series)

        else:
            frame.to_csv(summaryFileName, sep = '\t')
