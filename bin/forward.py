#!/usr/bin/env python3
import argparse
import csv
import sys
from Bio.Seq import Seq


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('seqtab', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    fieldnames = ['specimen', 'weight', 'seq']
    out = csv.DictWriter(args.out, fieldnames=fieldnames, lineterminator='\n')
    for r in csv.DictReader(args.seqtab, fieldnames=fieldnames):
        seq = Seq(r['seq'])
        r['seq'] = str(seq.reverse_complement())
        out.writerow(r)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
