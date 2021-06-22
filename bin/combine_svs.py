#!/usr/bin/env python3

"""Read vsearch output to determine SVs that are complements and
combine counts of complementary SVs on a per-sample basis

"""

import argparse
import csv


class DevNull:
    def write(*args, **kwargs):
        pass

    def writerow(*args, **kwargs):
        pass


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vsearch_out', help="output from vsearch task",
                        type=argparse.FileType('r'))
    parser.add_argument('weights', help="weights from write_seqs task",
                        type=argparse.FileType('r'))
    # output
    parser.add_argument('--corrected_weights', help="output combined weights of complement SVs per sample",
                        type=argparse.FileType('w'))


    args = parser.parse_args(arguments)
    vsearch_out = args.vsearch_out or DevNull()
    weights = args.weights or DevNull()
    corrected_weights = csv.writer(args.corrected_weights) if args.corrected_weights else DevNull()

    pass

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))