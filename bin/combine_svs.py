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
    parser.add_argument('weights', help="filename of weights output from write_seqs task",
                        type=str)
    # output
    parser.add_argument('--corrected_weights', help="filename for outputting corrected SV weights data, including combined complementary SVs",
                        type=str)


    args = parser.parse_args(arguments)
    vsearch_out = args.vsearch_out or DevNull()

    
    samples_per_label = defaultdict(list)
    count_per_sample = dict()

    with open(args.weights) as weights_file:
        weights = csv.DictReader(weights_file, fieldnames=['rep', 'sv', 'count'])

        for weight in weights:
            #weight['sample_id'] = weight['sv'].split(":")[1]
            samples_per_label[weight['rep']].append(weight['sv'])
            count_per_sample[weight['sv']] = weight['count']


    svs_with_complements = []
    corrected_counts = []

    complement_lines = [line.split() for line in vsearch_out]
    for line in complement_lines:
        # Use sv labels to infer which is more abundant and select that as 'primary' sv, other as complement to be combined
        if line[0].split(":")[0].lstrip("sv-") > line[1].split(":")[0].lstrip("sv-"):
            complement_sv = line[0]
            primary_sv = line[1]
            primary_sv_strand = "fwd"
            svs_with_complements.append(primary_sv)
            svs_with_complements.append(complement_sv)
        else:
            primary_sv = line[0]
            complement_sv = line[1]
            primary_sv_strand = "rev"
            svs_with_complements.append(primary_sv)
            svs_with_complements.append(complement_sv)

        complement_samples = samples_per_label[complement_sv]
        for complement_sample in complement_samples:
            complement_weight = count_per_sample[complement_sample]
            sample_id = complement_sample.split(":")[1]
            primary_samples = samples_per_label[primary_sv]
            for primary_sample in primary_samples:
                if primary_sample.split(":")[1] == sample_id:
                    summed_count = int(complement_weight) + int(count_per_sample[primary_sample])
                    corrected_counts.append({'rep':primary_sv, 'sv':primary_sample, 'count':summed_count, 'strand': primary_sv_strand})

    
    with open(args.weights) as weights_file:
        weights = csv.DictReader(weights_file, fieldnames=['rep', 'sv', 'count'])
        for weight in weights:
            if weight['rep'] not in svs_with_complements:
                corrected_counts.append({'rep': weight['rep'], 'sv': weight['sv'], 'count': weight['count'], 'strand':None})


    # TODO: get strand of svs from weights that weren't in vsearch output and carry them over to corrected_counts
    
    sorted_corrected_counts = sorted(corrected_counts, key=lambda k: int(k['count']), reverse=True)


    with open(args.corrected_weights, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames = ['rep', 'sv', 'count', 'strand'])
        writer.writerows(sorted_corrected_counts)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
