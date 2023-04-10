#!/usr/bin/env python3
import argparse
import csv
import fastalite
import sys


# Fixed as of version v2.18.0
def vsearch_issue_453(d):
    # strip off ;size=integer
    # https://github.com/torognes/vsearch/issues/453
    new_d = {}
    for k, v in d.items():
        k, v = k.split(';')[0], v.split(';')[0]
        new_d[k] = v
    return new_d


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clusters')
    parser.add_argument('fasta')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    with open(args.fasta) as fa_file:
        seeds = fastalite.fastalite(fa_file)
        seeds = (f.id.rsplit(';', 1) for f in seeds)
        seeds = {i: int(w.lstrip('size=')) for i, w in seeds}
    # merge Hit seqs' weights into the Seed seqs
    with open(args.clusters) as clusters_file:
        clusters = (row for row in clusters_file if row.startswith('H'))
        clusters = (row.split()[8:] for row in clusters)
        # clusters = vsearch_issue_453(clusters)
        for hit, seed in clusters:
            seeds[seed] += seeds[hit]
            del seeds[hit]  # remove non-seed cluster members
    out = csv.writer(args.out)
    with open(args.fasta) as fa_file:
        for f in fastalite.fastalite(fa_file):
            name = f.id.rsplit(';', 1)[0]
            if name in seeds:
                _, specimen = name.split(';')
                specimen = specimen.lstrip('specimen=')
                out.writerow([specimen, seeds[name], f.seq])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
