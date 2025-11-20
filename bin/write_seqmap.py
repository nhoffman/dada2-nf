#!/usr/bin/env python3

"""
Maps raw seqs to svs
"""

import sys
import argparse
import csv
import itertools


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sv_map')
    parser.add_argument('merged')
    parser.add_argument('r1')
    parser.add_argument('r2')
    parser.add_argument('out_merged')
    parser.add_argument('out_r1')
    parser.add_argument('out_r2')
    args = parser.parse_args(arguments)
    reader = csv.DictReader(
        open(args.sv_map),
        fieldnames=['specimen', 'direction', 'seqtab', 'representative']
        )
    svmap = {}
    for row in reader:
        key = (row['specimen'], row['direction'], row['seqtab'])
        svmap[key] = row['representative']
    directions = [['R1', args.r1, args.out_r1], ['R2', args.r2, args.out_r2]]
    for direction, rfn, ofn in directions:
        with open(ofn, 'w') as flo:
            out = csv.writer(flo)
            reader = csv.DictReader(
                open(rfn),
                fieldnames=['specimen', 'seqname', 'seqtab']
                )
            for row in reader:
                key = (row['specimen'], direction, row['seqtab'])
                sv = svmap.get(key, None)
                out.writerow([row['specimen'], row['seqname'], sv])
    with open(args.out_merged, 'w') as flo:
        out = csv.writer(flo)
        reader = csv.reader(open(args.merged))
        groupby = itertools.groupby(reader, key=lambda x: x[0])
        for specimen, grp in groupby:
            for i, (_, ir1, ir2) in enumerate(grp, 1):
                out.writerow([
                    svmap[(specimen, 'merged', str(i))],
                    svmap[(specimen, 'R1', ir1)],
                    svmap[(specimen, 'R2', ir2)]
                    ])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
