#!/usr/bin/env python3

"""Read vsearch output to determine SVs that are complements and
combine counts of complementary SVs on a per-sample basis

"""

import argparse
import csv
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

    lines = [line.split() for line in vsearch_out]

    weights = csv.reader(args.weights)
    counts = dict()
    for rep, sv, count in weights:
        counts[sv] = count

    for line in lines:
        # Assuming we should call the sv with lower numeric portion of sv name primary and other sv the complement
        if line[0].lstrip("sv-") < line[1].lstrip("sv-"):
            primary_sv = line[0]
            complement_sv = line[1]
        elif line[1].lstrip("sv-") < line[0].lstrip("sv-"):
            primary_sv = line[1]
            complement_sv = line[0]
        else:
            print("cannot determine primary vs complement sv from line in vsearch:")
            print(line)
            sys.exit(1)
        
        if primary_sv in counts:
            primary_count = counts[primary_sv]
            if complement_sv in counts:
                complement_count = counts[complement_sv]
                print("Primary sv: " + primary_sv)
                print("Primary sv count: " + str(primary_count))
                print("Complement sv: " + complement_sv)
                print("Complement sv count: " + str(complement_count))
                counts[primary_sv] = primary_count + complement_count
                counts.pop(complement_sv)

    print(len(counts))
    sorted_counts = {k:v for k, v in sorted(counts.items(), key=lambda x: x[1])}
    corrected_weights.writerow(["combinedSV", "count"])
    for k, v in sorted_counts.items():
        corrected_weights.writerow([k,v])
    


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
