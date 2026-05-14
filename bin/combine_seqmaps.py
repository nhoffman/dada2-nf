#!/usr/bin/env python3
import argparse
import csv
import sys


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clusters')
    parser.add_argument('fasta')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        help='three column seqtab [specimen,seed,sequence]',
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
