#!/usr/bin/env python3

"""
Required fields for the manifest:
* sampleid - a string found in the fastq file name uniquely identifying a
  specimen
* batch - a label grouping specimens into PCR batches for dada2::learnErrors()
* R1 - Path to the R1 FASTQ file
* R2 - Path to the R2 FASTQ file
* I1 - Path to the I1 index FASTQ file
* I2 - Path to the I2 index FASTQ file

Verifies the following:
* all sampleids are unique
* every sampleid in the manifest has corresponding R{1,2} and I{1,2}

"""

import argparse
import csv
import itertools
import operator
import sys

import openpyxl

KEEPCOLS = {'sampleid', 'sample_name', 'project', 'batch', 'controls',
            'R1', 'R2', 'I1', 'I2',
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


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--manifest', help="Manifest in excel or csv format")

    parser.add_argument('-b', '--batches',
                        help="Output csv file mapping sampleid to batch",
                        type=argparse.FileType('w'))
    parser.add_argument('-s', '--sample-info',
                        help="write the manifest as csv (requires -m/--manifest)",
                        type=argparse.FileType('w'))
    parser.add_argument('--sample-index',
                        help="Output csv file mapping index files to sampleid and R1/R2 files",
                        type=argparse.FileType('w'))
    parser.add_argument('--index-file-type', choices=['single', 'dual', 'none'],
                        default='dual', help='dual, single, or no index files?')
    args = parser.parse_args(arguments)

    # confirm that every sampleid is represented by the appropriate number of fastq files
    expected_labels = {
        'dual': ['R1', 'R2', 'I1', 'I2'],
        'single': ['R1', 'R2', 'I1'],
        'none': ['R1', 'R2'],
    }[args.index_file_type]

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
        assert len(manifest) == len(manifest_sampleids), "All sampleids must be unique"

        # confirm that all sampleids in the manifest have corresponding
        # fastq files
        extras = [
            row['sampleid']
            for row in manifest
            if any([
                row.get(cname) is None
                for cname in expected_labels
            ])
        ]
        if extras:
            sys.exit('samples in the manifest without fastq files: {}'.format(extras))

        if args.sample_info:
            writer = csv.DictWriter(args.sample_info, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(manifest)
    else:
        raise RuntimeError("Must provide --manifest")

    # Write out a CSV with just the sampleid and FASTQ paths in wide format
    if args.sample_index:
        out = csv.DictWriter(
            args.sample_index,
            fieldnames=['sampleid', 'direction', 'fastq', 'I1', 'I2'])
        for manifest_row in manifest:
            output_row = {'sampleid': manifest_row['sampleid']}

            if 'I1' in expected_labels:
                output_row['I1'] = manifest_row['I1']
            else:
                output_row['I1'] = ""
            if 'I2' in expected_labels:
                output_row['I2'] = manifest_row['I2']
            else:
                output_row['I2'] = ""

            for direction in ['R1', 'R2']:
                row['direction'] = direction
                row['fastq'] = manifest_row[direction]
                out.writerow(row)

    # finally, write an output file with columns (sampleid, batch)
    if args.batches:
        writer = csv.DictWriter(
            args.batches,
            fieldnames=['sampleid', 'batch'],
            extrasaction='ignore'
        )
        writer.writerows(manifest)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
