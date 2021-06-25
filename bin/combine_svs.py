#!/usr/bin/env python3

"""Read vsearch output to determine SVs that are complements and
combine counts of complementary SVs on a per-sample basis

"""

import argparse
from collections import defaultdict
import csv
import pandas as pd
import sys


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
    corrected_weights = csv.writer(args.corrected_weights) if args.corrected_weights else DevNull()


    weights = csv.DictReader(args.weights, fieldnames=['rep', 'sv', 'count'])
    
    samples_per_label = defaultdict(list)
    count_per_sample = dict()
    for weight in weights:
        weight['sample_id'] = weight['sv'].split(":")[1]
        samples_per_label[weight['rep']].append(weight['sv'])
        count_per_sample[weight['sv']] = weight['count']

    corrected_counts = []
    complement_lines = [line.split() for line in vsearch_out]
    for line in complement_lines:
        # Assume target SV (forward seq) is primary and query SV (reverse seq) is complement
        primary_sv = line[1]
        complement_sv = line[0]
        complement_samples = samples_per_label[complement_sv]
        for complement_sample in complement_samples:
            complement_weight = count_per_sample[complement_sample]
            sample_id = complement_sample.split(":")[1]
            primary_samples = samples_per_label[primary_sv]
            for primary_sample in primary_samples:
                if primary_sample.split(":")[1] == sample_id:
                    summed_count = int(complement_weight) + int(count_per_sample[primary_sample])
                    corrected_counts.append({'rep':primary_sv, 'sv':primary_sample, 'count':summed_count})

    print(corrected_counts)

    # TODO: get lines from weights that weren't in vsearch output and carry them over to corrected_counts
    # TODO: format corrected_counts properly into csv same shape as weights.csv input (rep, sv, count)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
