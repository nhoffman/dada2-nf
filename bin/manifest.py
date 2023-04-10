#!/usr/bin/env python3

"""
Required fields for the manifest:
* sampleid - a string found in the fastq file name uniquely identifying a
  specimen
* batch - a label grouping specimens into PCR batches for dada2::learnErrors()

Verifies the following:
* all sampleids are unique
* every sampleid in the manifest has corresponding R{1,2} and I{1,2}
* all fastq file names have sampleid as the first underscore-delimited token

"""

import argparse
import csv
import gzip
import itertools
import operator
import os
import sys
import re

import openpyxl

KEEPCOLS = {'sampleid', 'sample_name', 'project', 'batch', 'controls',
            'label', 'barcode_id', 'run_number', 'blast_database',
            'desc', 'other', 'quant', 'paired_ntc'}


def read_manifest_excel(fname, keepcols=KEEPCOLS):
    """Read the first worksheet from an excel file and return a generator
    of dicts with keys limited to 'keepcols'.

    """

    valgetter = operator.attrgetter('value')

    wb = openpyxl.load_workbook(fname)
    sheet = wb[wb.sheetnames[0]]

    rows = (r for r in sheet.iter_rows() if r[0].value)

    # get fieldnames from cells in the first row up to the first empty one
    header = itertools.takewhile(valgetter, next(rows))

    # replace column names to be discarded with '_'
    fieldnames = [(cell.value if cell.value in keepcols else '_')
                  for cell in header]
    popextra = '_' in fieldnames

    yield [f for f in fieldnames if f != '_']
    for row in rows:
        d = dict(zip(fieldnames, map(valgetter, row)))
        if popextra:
            d.pop('_')
        if d['sampleid']:
            yield d


def read_manifest_csv(fname, keepcols=KEEPCOLS):

    with open(fname) as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [(n if n in keepcols else '_')
                             for n in reader.fieldnames]

        yield [f for f in reader.fieldnames if f != '_']
        popextra = '_' in reader.fieldnames
        for d in reader:
            if popextra:
                d.pop('_')
            if d['sampleid']:
                yield dict(d)


def get_sampleid(pth):
    return os.path.basename(pth).split('_')[0]


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fastq_files', type=argparse.FileType('r'),
                        help="File listing fastq inputs")
    parser.add_argument('-m', '--manifest', help="Manifest in excel or csv format")

    parser.add_argument('-b', '--batches',
                        help="Output csv file mapping sampleid to batch",
                        type=argparse.FileType('w'))
    parser.add_argument('--counts',
                        help="raw fastq_file counts",
                        type=argparse.FileType('w'))
    parser.add_argument('-s', '--sample-info',
                        help="write the manifest as csv (requires -m/--manifest)",
                        type=argparse.FileType('w'))
    parser.add_argument('--index-file-type', choices=['single', 'dual', 'none'],
                        default='dual', help='dual, single, or no index files?')
    args = parser.parse_args(arguments)

    fq_files = sorted(line.strip() for line in args.fastq_files if line.strip())
    fq_sampleids = {get_sampleid(pth): pth for pth in fq_files}

    if args.manifest:
        if args.manifest.endswith('.csv'):
            read_manifest = read_manifest_csv
        else:
            read_manifest = read_manifest_excel

        manifest_reader = read_manifest(args.manifest)
        fieldnames = next(manifest_reader)
        manifest = list(manifest_reader)
        manifest_sampleids = {row['sampleid'] for row in manifest}

        # add missing batch labels
        for row in manifest:
            row['batch'] = row['batch'] or 'unknown'

        # make sure all sampleids are unique
        assert len(manifest) == len(manifest_sampleids)

        # confirm that all sampleids in the manifest have corresponding
        # fastq files
        extras = manifest_sampleids - fq_sampleids.keys()
        if extras:
            sys.exit('samples in the manifest without fastq files: {}'.format(extras))

        if args.sample_info:
            writer = csv.DictWriter(args.sample_info, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(manifest)
    else:
        manifest = [{'sampleid': sampleid, 'batch': 'unknown'}
                    for sampleid in sorted(fq_sampleids)]

    # confirm that every sampleid is represented by four fastq files
    expected_labels = {
        'dual': ['I1', 'I2', 'R1', 'R2'],
        'single': ['I1', 'R1', 'R2'],
        'none': ['R1', 'R2'],
    }[args.index_file_type]

    for sampleid, paths in itertools.groupby(fq_files, key=get_sampleid):
        labels = [re.findall(r'_([IR][12])_', fname)[0] for fname in paths]
        if labels != expected_labels:
            sys.exit('a fastq file is missing for sampleid {}: has {}'.format(
                sampleid, labels))

    # finally, write an output file with columns (sampleid, batch)
    if args.batches:
        writer = csv.DictWriter(
            args.batches, fieldnames=['sampleid', 'batch'], extrasaction='ignore')
        writer.writerows(manifest)

    if args.counts:
        writer = csv.DictWriter(args.counts, fieldnames=['sampleid', 'count'])
        writer.writeheader()
        for m in manifest:
            fq = gzip.open(os.path.basename(fq_sampleids[m['sampleid']]))
            count = sum(1 for li in fq if li.startswith(b'+'))
            writer.writerow({'sampleid': m['sampleid'], 'count': count})


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
