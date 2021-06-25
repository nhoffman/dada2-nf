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

            




    # weights = csv.reader(args.weights)
    # counts = dict()
    # for rep, sv, count in weights:
    #     label = rep + "," + sv
    #     counts[label] = count

    # weights = pd.read_csv(args.weights, header=None, names=["rep", "sv", "count"])
    # weights["sample_id"] = weights.apply(lambda row: (row.sv).split(":")[1], axis=1)

    # complement_lines = [line.split() for line in vsearch_out]
    # for line in complement_lines:
    #     # Assume target SV (forward seq) is primary and query SV (reverse seq) is complement
    #     primary_sv = line[1]
    #     complement_sv = line[0]
    #     primary_sv_weights = weights.loc[weights['rep'] == primary_sv]
    #     complement_sv_weights = weights.loc[weights['rep'] == complement_sv]
    #     import pdb;pdb.set_trace()

    # complements = dict()
    # for line in complement_lines:
    #     # Assume target SV (forward seq) is primary and query SV (reverse seq) is complement
    #     primary_sv = line[1]
    #     complement_sv = line[0]
    #     complements[primary_sv] = complement_sv
    
    # for primary_sv, complement_sv in complements.items():
    #     primary_sv_weights = weights.loc[weights['rep'] == primary_sv]
    #     complement_sv_weights = weights.loc[weights['rep'] == complement_sv]
    #     import pdb;pdb.set_trace()

    # corrected_weights = dict()

    # weights = csv.DictReader(args.weights, fieldnames=['rep', 'sv', 'count'])
    # for row in weights:
    #     if row['rep'] in complements.keys() 
    
        
    #     if primary_sv in counts:
    #         primary_count = counts[primary_sv]
    #         if complement_sv in counts:
    #             complement_count = counts[complement_sv]
    #             print("Primary sv: " + primary_sv)
    #             print("Primary sv count: " + str(primary_count))
    #             print("Complement sv: " + complement_sv)
    #             print("Complement sv count: " + str(complement_count))
    #             counts[primary_sv] = primary_count + complement_count
    #             counts.pop(complement_sv)

    # print(len(counts))
    # sorted_counts = {k:v for k, v in sorted(counts.items(), key=lambda x: x[1])}
    # corrected_weights.writerow(["combinedSV", "count"])
    # for k, v in sorted_counts.items():
    #     corrected_weights.writerow([k,v])
    


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
