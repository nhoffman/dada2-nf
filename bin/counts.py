#!/usr/bin/env python3
import argparse
import csv
import itertools
import sys

STEP_ORDER = [
    'raw',
    'barcodecop',
    'cutadapt',
    'split',
    'dada2_filtered_and_trimmed',
    'dada2_denoised',
    'dada2']


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', type=argparse.FileType('r'))
    parser.add_argument('cutadapt', type=argparse.FileType('r'))
    parser.add_argument('split_orientations', type=argparse.FileType('r'))
    parser.add_argument('barcodecop', type=argparse.FileType('r'))
    parser.add_argument('dada2', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    rows = []
    for m in csv.DictReader(args.manifest):
        m['step'] = 'raw'
        rows.append(m)
    barcodecop = csv.DictReader(
        args.barcodecop,
        fieldnames=['sampleid', 'in', 'count'])
    bcop = {}
    for b in barcodecop:
        # deduplicate bcop R1/R2 rows
        sampleid = b['sampleid'].split('_')[0]
        b['sampleid'] = sampleid
        b['step'] = 'barcodecop'
        bcop[sampleid] = b
    rows.extend(bcop.values())
    for c in csv.DictReader(args.cutadapt):
        c['count'] = c['out_reads']
        c['step'] = 'cutadapt'
        rows.append(c)
    splits = csv.DictReader(
        args.split_orientations,
        fieldnames=['sampleid', 'orientation', 'count'])
    for s in splits:
        s['step'] = 'split'
        rows.append(s)
    dada2 = list(csv.DictReader(args.dada2))
    for d in dada2:
        d['step'] = 'dada2_filtered_and_trimmed'
        d['count'] = d['filtered_and_trimmed']
        rows.append(d.copy())
    for d in dada2:
        d['step'] = 'dada2_denoised'
        d['direction'] = 'R1'
        d['count'] = d['denoised_r1']
        rows.append(d.copy())
        d['step'] = 'dada2_denoised'
        d['direction'] = 'R2'
        d['count'] = d['denoised_r2']
        rows.append(d.copy())
    for d in dada2:
        d['step'] = 'dada2'
        d['direction'] = 'merged'
        d['count'] = d['merged']
        rows.append(d.copy())
    out = csv.DictWriter(
        args.out,
        fieldnames=['step', 'sampleid', 'orientation', 'direction', 'count', 'total'],
        extrasaction='ignore')
    rows = sorted(rows, key=lambda x: (x['step'], x['sampleid']))
    counts = {}
    for k, g in itertools.groupby(rows, key=lambda x: (x['step'], x['sampleid'])):
        counts[k] = sum(int(r['count']) for r in g)
    for r in rows:
        r['total'] = counts[(r['step'], r['sampleid'])]
    rows = sorted(rows, key=lambda x: (int(x['total']), int(x['count'])), reverse=True)
    rows = sorted(rows, key=lambda x: STEP_ORDER.index(x['step']))
    out.writeheader()
    out.writerows(rows)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
