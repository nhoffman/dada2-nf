#!/usr/bin/env python3
"""
Parse miseq fastqs into a dada2-nf compatible manifest

Output csv:
    sampleid,project,batch,datadir,R1,R2,I1,I2
"""
import argparse
import csv
import glob
import itertools
import os
import re
import sys


# FASTQ sample IDs are parsed in two steps:
# 1. READ_SUFFIX removes the read/index marker and everything after it.
# 2. SAMPLE_SUFFIX removes optional Illumina sample/lane suffixes.
#
# Examples:
#   m3n701-s502_S1_L001_R1_001.fastq.gz -> m3n701-s502
#   22R255-NGSITS49_S4_L001_I2_001.fastq.gz -> 22R255-NGSITS49
#   1_07_5479_wk12_R2_001.fastq.gz -> 1_07_5479_wk12
SAMPLE_SUFFIX = re.compile(r'_S\d+(_L\d{3})?$')
READ_SUFFIX = re.compile(r'_[RI][12]_.*$')


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fastqs', help='Fastq directory or fastq_list.txt')
    parser.add_argument(
        '--project',
        help='Value for optional project column')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    if args.project:
        fieldnames = ['sampleid', 'project', 'batch',
                      'datadir', 'R1', 'R2', 'I1', 'I2']
    else:
        fieldnames = ['sampleid', 'batch',
                      'datadir', 'R1', 'R2', 'I1', 'I2']
    out = csv.DictWriter(
        args.out,
        extrasaction='ignore',
        fieldnames=fieldnames
        )
    out.writeheader()
    fastqs = read_fastqs(args.fastqs)
    fastqs = itertools.groupby(sorted(fastqs, key=samplename), key=samplename)
    for sampleid, fls in fastqs:
        sampleid = SAMPLE_SUFFIX.sub('', sampleid)
        fls = list(fls)
        datadir = os.path.commonpath([os.path.dirname(f) for f in fls])
        reads = read_files(fls, datadir)
        out.writerow({
            'sampleid': sampleid,
            'project': args.project,
            'batch': 1,
            'datadir': os.path.abspath(datadir),
            'R1': reads['R1'],
            'R2': reads['R2'],
            'I1': reads['I1'],
            'I2': reads['I2']
            })


def read_fastqs(path):
    if os.path.isdir(path):
        return glob.glob(os.path.join(path, '*.fastq.gz'))
    with open(path) as handle:
        return [
            line.strip() for line in handle
            if line.strip().endswith('.fastq.gz')
        ]


def read_files(fastqs, datadir):
    reads = {'R1': None, 'R2': None, 'I1': None, 'I2': None}
    for fastq in fastqs:
        basename = os.path.basename(fastq)
        for read in reads:
            if '_' + read + '_' in basename:
                reads[read] = os.path.relpath(fastq, datadir)
                break
        else:
            raise ValueError('unknown read in filename: ' + fastq)
    return reads


def samplename(fl):
    basename = os.path.basename(fl)
    sample = READ_SUFFIX.sub('', basename)
    if sample == basename:
        raise ValueError('unknown sample filename: ' + fl)
    return sample


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
