#!/usr/bin/env python3

"""Read alignment scores and indicate sequences with a score
above a specified score.
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
        strand = colnames.index('strand')
        aligns = (r.split() for r in args.scores if not r.startswith('#'))
        aligns = {r[name]: r for r in aligns}
    elif args.strategy == "vsearch":
        colnames = ['query', 'id', 'qstrand']
        name = colnames.index('query')
        score = colnames.index('id')
        strand = colnames.index('qstrand')
        aligns = (r.split() for r in args.scores)
        aligns = {r[name]: r for r in aligns}
    else:
        aligns = {}

    def forward(seq):
        if seq.id in aligns and aligns[seq.id][strand] == '-':
            seq.seq = seq.seq.reverse_complement()
        return seq

    seqs = (forward(s) for s in SeqIO.parse(args.seqs, format='fasta'))

    outcomes.writerow(['seqname', 'is_target'])

    for s in seqs:
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
        f, r, n = defaultdict(int), defaultdict(int), defaultdict(int)
        for rep, sv, count in weights:
            specimen = sv.split(':')[-1]
            if args.strategy == 'none':
                f[specimen] += int(count)
            elif rep not in aligns:
                n[specimen] += int(count)
            elif float(aligns[rep][score]) < args.min_score:
                n[specimen] += int(count)
            elif aligns[rep][strand] == '+':
                f[specimen] += int(count)
            elif aligns[rep][strand] == '-':
                r[specimen] += int(count)
            else:
                err = 'unknown strand char: ' + aligns[rep][strand]
                raise ValueError(err)

        writer = csv.writer(args.counts)
        writer.writerow(['sampleid', 'target_f', 'target_r', 'not_target'])
        for specimen in sorted(set(f) | set(r) | set(n)):
            writer.writerow([specimen, f[specimen], r[specimen], n[specimen]])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
