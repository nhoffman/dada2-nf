#!/usr/bin/env python3
"""
Output sequences into forward, reverse and off-target
"""
import argparse
import csv
import gzip
import os
import sys

from Bio import SeqIO


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sampleid')
    parser.add_argument('orientations')
    parser.add_argument('sample', nargs='+')
    parser.add_argument('--counts')
    parser.add_argument('--fwd-out', default='forward')
    parser.add_argument('--rev-out', default='reverse')
    parser.add_argument('--off-out', default='off_target')
    args = parser.parse_args(arguments)
    os.makedirs(args.fwd_out, exist_ok=True)
    os.makedirs(args.rev_out, exist_ok=True)
    os.makedirs(args.off_out, exist_ok=True)
    orientations = dict(csv.reader(open(args.orientations), delimiter='\t'))
    for reads in args.sample:
        name = os.path.basename(reads)
        fwd_count = 0
        rev_count = 0
        off_count = 0
        with (gzip.open(os.path.join(args.fwd_out, name), 'wt') as fout,
              gzip.open(os.path.join(args.rev_out, name), 'wt') as rout,
              gzip.open(os.path.join(args.off_out, name), 'wt') as oout):
            for r in SeqIO.parse(gzip.open(reads, 'rt'), 'fastq'):
                if r.name not in orientations:
                    SeqIO.write(r, oout, 'fastq')
                    off_count += 1
                elif orientations[r.name] == '+':
                    SeqIO.write(r, fout, 'fastq')
                    fwd_count += 1
                elif orientations[r.name] == '-':
                    SeqIO.write(r, rout, 'fastq')
                    rev_count += 1
                else:
                    raise ValueError(
                        'Unknown orientation: '
                        f'{r.name} {orientations[r.name]}')
    if args.counts:
        with open(args.counts, 'w') as out:
            counts_out = csv.DictWriter(
                out,
                fieldnames=['sampleid', 'orientation', 'reoriented'])
            counts_out.writerow({
                'sampleid': args.sampleid,
                'orientation': 'forward',
                'reoriented': fwd_count})
            counts_out.writerow({
                'sampleid': args.sampleid,
                'orientation': 'reverse',
                'reoriented': rev_count})
            counts_out.writerow({
                'sampleid': args.sampleid,
                'orientation': 'off_target',
                'reoriented': off_count})


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
