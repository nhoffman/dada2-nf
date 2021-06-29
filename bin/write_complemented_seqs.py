#!/usr/bin/env python3


import argparse
from Bio.Seq import Seq
import csv
from fastalite import fastalite
import os
import sys


class DevNull:
    def write(self, *args, **kwargs):
        pass

    def writerow(self, *args, **kwargs):
        pass


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'seqs', type=argparse.FileType('r'),
        help='seqs.fasta file containing all seqs, output from initial write_seqs task')
    parser.add_argument(
        'corrected_weights', type=str,
        help='corrected_weights.csv file output from combine_svs task')
    parser.add_argument(
        '--final_seqs', type=argparse.FileType('w'),
        help=('fasta file for outputting the final sequences, with complements removed ',
        'and all seqs in forward direction'), default='final_complemented_seqs.fasta')

    args = parser.parse_args(arguments)


    final_seqs = args.final_seqs or DevNull()

    orig_seqs = dict()
    orig_sv_labels = []
    for orig_seq in fastalite(args.seqs):
        orig_seqs[orig_seq.id] = orig_seq.seq
        orig_sv_labels.append(orig_seq.id)


    final_sv_labels = []
    svs_in_weights_not_in_orig = []
    with open(args.corrected_weights) as weights_file:
        weights = csv.DictReader(weights_file, fieldnames=['rep', 'sv', 'count', 'strand', 'merged'])
        for weight in weights:
            label = weight['sv']
            final_sv_labels.append(label)
            if weight['strand'] == 'fwd':
                sv_seq = orig_seqs.get(label)
                if sv_seq is not None:
                    final_seqs.write('>{}\n{}\n'.format(label, sv_seq))
                else:
                    svs_in_weights_not_in_orig.append(label)
            elif weight['strand'] == 'rev':
                sv_seq = orig_seqs.get(label)
                if sv_seq is not None:
                    rev_comp = str(Seq(sv_seq).reverse_complement())
                    final_seqs.write('>{}\n{}\n'.format(label, rev_comp))
                else:
                    svs_in_weights_not_in_orig.append(label)

    print("final_sv_labels")
    print(len(final_sv_labels))
    print("orig_sv_labels")
    print(len(orig_seqs))
    print("svs in weight not in orig")
    print(svs_in_weights_not_in_orig)
    print(len(svs_in_weights_not_in_orig))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

