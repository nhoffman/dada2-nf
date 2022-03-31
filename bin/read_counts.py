#!/usr/bin/env python3
import argparse
import csv
import itertools
import sys

from fastalite import fastqlite, Opener


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'fastq',
        help='reads to count in fastq format',
        metavar='file.fastq[.bz2|.gz]',
        type=Opener(),
        )
    parser.add_argument(
        'fastq_out',
        help='output fastq',
        type=Opener('w'),
        )
    parser.add_argument(
        'read_counts',
        help='tabulate read counts and store as a CSV',
        metavar='FILE',
        type=argparse.FileType('w'),
        )
    parser.add_argument(
        '--head',
        help='limit the output file to N records',
        metavar='N',
        type=int)
    args = parser.parse_args(arguments)
    fastq = itertools.islice(fastqlite(args.fastq), args.head)
    count = 0
    for count, f in enumerate(fastq, start=1):
        args.fastq_out.write(f'@{f.description}\n{f.seq}\n+\n{f.qual}\n')
    read_counts_writer = csv.writer(args.read_counts)
    read_counts_writer.writerow([args.fastq_out.name, count, count])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
