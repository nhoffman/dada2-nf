#!/usr/bin/env python3
"""
Parse miseq fastqs into a dada2-nf compatiable manifest

Output csv:
    sampleid,project,batch,datadir,R1,R2,I1,I2
"""
import argparse
import csv
import glob
import itertools
import os
import sys


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('datadir', help='Fastq directory')
    parser.add_argument(
        '--project',
        help='Value for optional project column')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    if args.project:
        fieldnames = ['sampleid', 'project', 'batch',
                      'datadir', 'R1', 'R2', 'I1', 'I2']
    else:
        fieldnames = ['sampleid', 'batch',
                      'datadir', 'R1', 'R2', 'I1', 'I2']
    out = csv.DictWriter(
        args.out,
        extrasaction='ignore',
        fieldnames=fieldnames
        )
    out.writeheader()
    fastqs = glob.glob(os.path.join(args.datadir, '*.fastq.gz'))
    fastqs = (os.path.basename(f) for f in fastqs)
    fastqs = itertools.groupby(sorted(fastqs), key=samplename)
    for sampleid, fls in fastqs:
        # rename sampleid with the groupby _S
        sampleid = sampleid.split('_')[0]
        R1, R2, *index = list(fls)
        if len(index) == 2:
            I1, I2 = index
        elif len(index) == 1:
            I1, = index
            I2 = None
        else:
            I1, I2 = None, None
        out.writerow({
            'sampleid': sampleid,
            'project': args.project,
            'batch': 1,
            'datadir': os.path.abspath(args.datadir),
            'R1': R1,
            'R2': R2,
            'I1': I1,
            'I2': I2
            })


def samplename(fl):
    if '_R' in fl:
        return os.path.basename(fl).split('_R')[0]
    elif '_I' in fl:
        return os.path.basename(fl).split('_I')[0]
    else:
        raise ValueError('unknown sample filename: ' + fl)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
