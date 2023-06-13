#!/usr/bin/env python3
import argparse
import csv
import os
import sys


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    for row in csv.DictReader(args.manifest):
        datadir = row['datadir'].strip()
        for col in ['R1', 'R2', 'I1', 'I2']:
            fq = row[col].strip()
            if fq:
                args.out.write(os.path.join(datadir, fq) + '\n')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
