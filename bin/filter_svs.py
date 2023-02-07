#!/usr/bin/env python3

"""Read alignment scores and indicate sequences with a score
above a specified score.

NOTE: sequences filtered out by `cmsearch -E <x>` or `vsearch` and
not included in the alignment scores input will still be reported
in the `--failing` output
"""

import sys
import argparse
import csv
from collections import defaultdict

from Bio import SeqIO


class DevNull:
    def write(*args, **kwargs):
        pass

    def writerow(*args, **kwargs):
        pass


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('seqs', help="input sequences",
                        type=argparse.FileType('r'))
    parser.add_argument('orientation', help="original sequence orientation")
    parser.add_argument('scores', help="strategy alignment scores",
                        type=argparse.FileType('r'))
    parser.add_argument('--weights',
                        type=argparse.FileType('r'))
    # outputs
    parser.add_argument('--passing', help="output fasta of passing sequences",
                        type=argparse.FileType('w'))
    parser.add_argument('--failing', help="outpt fasta of failing sequences",
                        type=argparse.FileType('w'))
    parser.add_argument('--outcomes', help="output outcome for each sv",
                        type=argparse.FileType('w'))
    parser.add_argument('--counts', type=argparse.FileType('w'),
                        help=("output counts of 'target' and 'other' reads "
                              "per specimen (requires --weights)"))
    parser.add_argument('--min-score', type=int, default=0,
                        help='minimum bit score [default %(default)s]')
    parser.add_argument('--strategy',
                        choices=['cmsearch', 'vsearch', 'none'],
                        default='cmsearch',
                        help=("script strategy for alignments "
                              "scores and orientation [default %(default)s]"))

    args = parser.parse_args(arguments)
    passing = args.passing or DevNull()
    failing = args.failing or DevNull()
    outcomes = csv.writer(args.outcomes) if args.outcomes else DevNull()

    if args.strategy == "cmsearch":
        colnames = ['target name', 'target accession', 'query name',
                    'query accession', 'mdl', 'mdl from', 'mdl to',
                    'seq from', 'seq to', 'strand', 'trunc', 'pass',
                    'gc', 'bias', 'score', 'e_value', 'inc',
                    'description of target']
        name = colnames.index('target name')
        score = colnames.index('score')
        aligns = (r.split() for r in args.scores if not r.startswith('#'))
        aligns = {r[name]: r for r in aligns}
    elif args.strategy == "vsearch":
        colnames = ['query', 'id', 'qstrand']
        name = colnames.index('query')
        score = colnames.index('id')
        aligns = (r.split() for r in args.scores)
        aligns = {r[name]: r for r in aligns}
    else:
        aligns = {}

    outcomes.writerow(['seqname', 'is_target'])

    for s in SeqIO.parse(args.seqs, format='fasta'):
        output = '>{seq.id}\n{seq.seq}\n'.format(seq=s)
        if args.strategy == 'none':
            outcomes.writerow([s.id, True])
            passing.write(output)
        elif s.id in aligns and float(aligns[s.id][score]) >= args.min_score:
            outcomes.writerow([s.id, True])
            passing.write(output)
        else:
            outcomes.writerow([s.id, False])
            failing.write(output)

    if args.counts and args.weights:
        weights = csv.reader(args.weights)
        p, f = defaultdict(int), defaultdict(int)
        for rep, sv, count in weights:
            specimen = sv.split(':')[-1]
            if args.strategy == 'none':
                p[specimen] += int(count)
            elif rep not in aligns:
                f[specimen] += int(count)
            elif float(aligns[rep][score]) < args.min_score:
                f[specimen] += int(count)
            else:
                p[specimen] += int(count)

        writer = csv.writer(args.counts)
        writer.writerow(['sampleid', 'orientation', 'target', 'not_target'])
        for specimen in sorted(set(p) | set(f)):
            row = [specimen, args.orientation, p[specimen], f[specimen]]
            writer.writerow(row)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
