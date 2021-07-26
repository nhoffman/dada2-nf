#!/usr/bin/env python3
"""
Appends size/weight/abundance information to
sequence header for vsearch --cluster_size
"""
import argparse
import collections
import csv
import fastalite
import sys


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('table', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    abundances = collections.defaultdict(int)
    for row in csv.DictReader(args.table):
        abundances[row['representative']] += int(row['count'])
    for f in fastalite.fastalite(args.fasta):
        header = '{};size={}'.format(f.id, abundances[f.id])
        args.out.write('>{}\n{}\n'.format(header, f.seq))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
