#!/usr/bin/env python3

"""Write counts of sequence variants.

Description of outputs:

# specimen_map.csv
sv-0001:m76n710-s511  m76n710-s511
sv-0001:m76n712-s506  m76n712-s506
sv-0001:m76n712-s505  m76n712-s505
sv-0001:m76n712-s511  m76n712-s511

# weights.csv
sv-0001:m76n710-s511  sv-0001:m76n710-s511  194200
sv-0001:m76n710-s511  sv-0001:m76n712-s506  169784
sv-0001:m76n710-s511  sv-0001:m76n712-s505  124221
sv-0001:m76n710-s511  sv-0001:m76n712-s511  110659

# dada2_sv_table.csv
sv                    m76n701-s502  m76n701-s503  m76n701-s505
sv-0001:m76n710-s511  27            94            122
sv-0002:...           36            31704         8829
sv-0003:...           0             8             0
sv-0004:...           0             0             34

# dada2_sv_table_long.csv
specimen      count   sv       representative
m76n710-s511  194200  sv-0001  sv-0001:m76n710-s511
m76n712-s506  169784  sv-0001  sv-0001:m76n710-s511
m76n712-s505  124221  sv-0001  sv-0001:m76n710-s511
m76n712-s511  110659  sv-0001  sv-0001:m76n710-s511

"""

import os
import sys
import argparse


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
        'reverse_seqs', type=argparse.FileType('w'),
        help='fasta file containing reads determined to be from reverse strand')
    parser.add_argument(
        'corrected_weights', type=argparse.FileType('r'),
        help='corrected_weights.csv file output from combine_svs task')
    parser.add_argument(
        '--final_seqs', type=argparse.FileType('w'),
        help='fasta file for outputting the final sequences, with complements removed ',
        'and all seqs in forward direction', default='final_complemented_seqs.fasta')

    args = parser.parse_args(arguments)


    orig_seqs = args.seqs or DevNull()
    rev_seqs = args.reverse_seqs or DevNull()
    weights = args.corrected_weights or DevNull()

    
    # TODO: reverse complement all seqs from rev_seqs so that all final_seqs are in forward direction
    # TODO: drop all seqs from final_seqs that are not in weights input, assuming they've been combined with their complement in previous task


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

