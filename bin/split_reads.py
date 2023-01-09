#!/usr/bin/env python3
"""
Split reads according to list of seqnames in file or fastq
"""
import argparse
import gzip
import os
import sys

from Bio import SeqIO


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sync', nargs='+')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--seqname-file')
    inputs.add_argument('--fastq')
    parser.add_argument('--outdir', default='passed')
    parser.add_argument('--out_suffix', default='')
    parser.add_argument('--dropdir', default='dropped')
    parser.add_argument('--drop_suffix', default='')
    args = parser.parse_args(arguments)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.dropdir, exist_ok=True)
    if args.seqname_file:
        keep = (s.strip() for s in open(args.seqname_file))
        keep = {k for k in keep if k}
    else:
        keep = {f.name for f in SeqIO.parse(args.fastq, 'fastq')}
    for reads in args.sync:
        name, ext = os.path.basename(reads).split('.', 1)
        outname = f'{name}{args.out_suffix}.{ext}'
        out = os.path.join(args.outdir, outname)
        dropname = f'{name}{args.drop_suffix}.{ext}'
        dropped = os.path.join(args.dropdir, dropname)
        with (gzip.open(out, 'wt') as out,
                gzip.open(dropped, 'wt') as dropped):
            for r in SeqIO.parse(gzip.open(reads, 'rt'), 'fastq'):
                if r.name in keep:
                    SeqIO.write(r, out, 'fastq')
                else:
                    SeqIO.write(r, dropped, 'fastq')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
