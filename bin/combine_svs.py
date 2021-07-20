#!/usr/bin/env python3
import argparse
import collections
import csv
import fastalite
import sys


def vsearch_issue_453(d):
    # strip off ;size=integer
    # https://github.com/torognes/vsearch/issues/453
    for k, v in list(d.items()):
        k, v = k.split(';')[0], v.split(';')[0]
        d[k] = v
    return d


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('table', type=argparse.FileType('r'))
    parser.add_argument('clusters', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    clusters = (row for row in args.clusters if row.startswith('H'))
    complements = dict(row.split()[8:] for row in clusters)
    complements = vsearch_issue_453(complements)
    reps = collections.defaultdict(lambda: collections.defaultdict(int))
    for row in csv.DictReader(args.table):
        rep = complements.get(row['representative'], row['representative'])
        reps[rep][row['specimen']] += int(row['count'])
    for f in fastalite.fastalite(args.fasta):
        for specimen, count in reps[f.id].items():
            args.out.write('{},{},{}\n'.format(specimen, count, f.seq))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
