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

import sys
import argparse
import csv
import math
from Bio.Seq import Seq
from itertools import groupby
from collections import defaultdict, OrderedDict
from operator import itemgetter
from multiprocessing import Pool


def read_seqtab(fname):
    with open(fname) as f:
        reader = csv.reader(f)
        return [(specimen, int(count), seq) for specimen, count, seq in reader]


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
        'seqtabs', nargs='*',
        help='One or more headerless CSV files with columns specimen,count,seq')
    parser.add_argument(
        '--seqtablist', type=argparse.FileType('r'),
        help='A file listing one seqtab file per line')

    parser.add_argument(
        '--seqs', type=argparse.FileType('w'),
        help='fasta file containing sequence variants')
    parser.add_argument(
        '--weights', type=argparse.FileType('w'),
        help='csv file with columns sv,sv-specimen,weight')
    parser.add_argument(
        '--specimen-map', type=argparse.FileType('w'),
        help='csv file with columns seqname,specimen')
    parser.add_argument(
        '--sv-table', type=argparse.FileType('w'),
        help='csv file with svs in rows and specimens in columns')
    parser.add_argument(
        '--sv-table-long', type=argparse.FileType('w'),
        help=('"long" format csv file with columns '
              'specimen,count,sv,representative'))
    parser.add_argument(
        '--specimen-table', type=argparse.FileType('w'),
        help='specimen counts')
    parser.add_argument(
        '--direction',
        help='label to add to all sequence names')
    parser.add_argument(
        '-j', '--num-processes', type=int, default=1,
        help="number of processes for reading files [%(default)s]")
    parser.add_argument(
        '--reverse-complement', action='store_true',
        help='reverse complement SVs')

    args = parser.parse_args(arguments)

    if args.seqtabs:
        seqtabfiles = args.seqtabs
    elif args.seqtablist:
        seqtabfiles = [line.strip() for line in args.seqtablist if line.strip()]
    else:
        sys.exit('Error: must specify one or more seqtab files')

    with Pool(processes=args.num_processes) as pool:
        rows = [r for f in pool.map(read_seqtab, seqtabfiles) for r in f]

    seqfile = args.seqs or DevNull()
    specimen_map = csv.writer(args.specimen_map) if args.specimen_map else DevNull()
    weights = csv.writer(args.weights) if args.weights else DevNull()
    sv_table = csv.writer(args.sv_table) if args.sv_table else DevNull()
    sv_table_long = csv.writer(args.sv_table_long) if args.sv_table_long else DevNull()

    svlist = []
    all_specimens = set()
    for seq, grp in groupby(sorted(rows, key=itemgetter(2, 1)), itemgetter(2)):
        # each group is ordered by count desc
        specimens, counts, __ = zip(*reversed(list(grp)))
        svlist.append((sum(counts), OrderedDict(zip(specimens, counts)), seq))
        all_specimens |= set(specimens)

    all_specimens = sorted(all_specimens)

    # order by overall count, desc; include hash of seq to make sure
    # sorting of ties is stable
    svlist.sort(key=lambda r: (r[0], hash(r[2])), reverse=True)
    padchars = math.ceil(math.log10(len(svlist) + 1))

    def svname(i, specimen=None):
        name = ['sv-{:0{}}'.format(i, padchars)]
        if specimen:
            name.append(specimen)
        return ':'.join(name)

    sv_table.writerow(['sv'] + all_specimens)
    sv_table_long.writerow(['specimen', 'count', 'sv', 'representative'])

    for i, (total, specimens, seq) in enumerate(svlist, 1):
        first_specimen = next(iter(specimens.keys()))
        representative = svname(i, first_specimen)

        if args.reverse_complement:
            seq = str(Seq(seq).reverse_complement())

        seqfile.write('>{}\n{}\n'.format(representative, seq))
        sv_table.writerow([representative] + \
                          [specimens.get(s, 0) for s in all_specimens])

        for specimen, count in specimens.items():
            this_seqname = svname(i, specimen)
            specimen_map.writerow([this_seqname, specimen])
            weights.writerow([representative, this_seqname, count])
            sv_table_long.writerow([specimen, count, svname(i), representative])

    # specimen counts
    counts = defaultdict(int)
    for specimen, count, _ in rows:
        counts[specimen] += count
    out = csv.DictWriter(
        args.specimen_table,
        fieldnames=['sampleid', 'direction', 'count'])
    for k, v in counts.items():
        out.writerow({'sampleid': k, 'direction': args.direction, 'count': v})


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
