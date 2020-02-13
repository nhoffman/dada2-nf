#!/usr/bin/env python3
import argparse
import csv
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
        'read_counts',
        help='tabulate read counts and store as a CSV',
        metavar='FILE',
        type=argparse.FileType('w'),
        )
    args = parser.parse_args(arguments)
    count = sum(1 for _ in fastqlite(args.fastq))
    read_counts_writer = csv.writer(args.read_counts)
    read_counts_writer.writerow([args.fastq.name, count, count])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
