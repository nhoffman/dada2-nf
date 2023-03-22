#!/usr/bin/env python3
import argparse
import csv
import sys


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('seqtabs', nargs='+', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    fieldnames = ['specimen', 'weight', 'seq']
    seqtabs = (csv.DictReader(s, fieldnames=fieldnames) for s in args.seqtabs)
    seqtabs = (r for s in seqtabs for r in s)  # flatten
    for i, s in enumerate(seqtabs):
        header = '{};specimen={};size={}'.format(i, s['specimen'], s['weight'])
        args.out.write('>{}\n{}\n'.format(header, s['seq']))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
