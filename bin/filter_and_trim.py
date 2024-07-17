#!/usr/bin/env python3

import argparse
import fastalite
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sampleid')
    parser.add_argument('orientation')

    parser.add_argument('r1', type=fastalite.Opener())
    parser.add_argument('r2', type=fastalite.Opener())

    parser.add_argument('r1_passed', type=fastalite.Opener())
    parser.add_argument('r2_passed', type=fastalite.Opener())

    outs = parser.add_argument_group(title='outputs')
    outs.add_argument(
        'r1_dropped', type=fastalite.Opener('w'), default=sys.stdout)
    outs.add_argument(
        'r2_dropped', type=fastalite.Opener('w'), default=sys.stdout)

    outs.add_argument(
        'counts', type=fastalite.Opener('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    passed = fastalite.fastqlite(args.r1_passed, allow_empty=True)
    passed = set(f.id for f in passed)
    for fq, out in [(args.r1, args.r1_dropped), (args.r2, args.r2_dropped)]:
        dropped = fastalite.fastqlite(fq, allow_empty=True)
        dropped = (f for f in dropped if f.id not in passed)
        for f in dropped:
            out.write(f'@{f.id}\n{f.seq}\n+\n{f.qual}\n')
    args.counts.write(
        'sampleid,orientation,count\n'
        '{},{},{}\n'.format(
            args.sampleid, args.orientation, len(passed))
        )


if __name__ == '__main__':
    sys.exit(main())
