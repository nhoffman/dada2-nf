#!/usr/bin/env python3
import argparse
import csv
import sys

STEP_ORDER = [
    'raw',
    'barcodecop',
    'cutadapt',
    'split',
    'filter_and_trim',
    'dada2',
    'svs']


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', type=argparse.FileType('r'))
    parser.add_argument('cutadapt', type=argparse.FileType('r'))
    parser.add_argument('split_orientations', type=argparse.FileType('r'))
    parser.add_argument('barcodecop', type=argparse.FileType('r'))
    parser.add_argument('dada2', type=argparse.FileType('r'))
    parser.add_argument('specimens', type=argparse.FileType('r'))
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
        d['step'] = 'filter_and_trim'
        d['count'] = d['filtered_and_trimmed']
        rows.append(d.copy())
    for d in dada2:
        d['step'] = 'dada2'
        d['direction'] = 'R1'
        d['count'] = d['denoised_r1']
        rows.append(d.copy())
        d['direction'] = 'R2'
        d['count'] = d['denoised_r2']
        rows.append(d.copy())
        d['direction'] = 'merged'
        d['count'] = d['no_chimeras']
        rows.append(d.copy())
    fieldnames = ['sampleid', 'direction', 'count']
    for s in csv.DictReader(args.specimens, fieldnames=fieldnames):
        rows.append({'step': 'svs', **s})
    fieldnames = ['step', 'sampleid', 'orientation', 'direction', 'count']
    for r in rows:
        if 'direction' not in r:
            r['direction'] = ''
        if 'orientation' not in r:
            r['orientation'] = ''
    out = csv.DictWriter(
        args.out,
        fieldnames=fieldnames,
        extrasaction='ignore')
    rows = sorted(
        rows,
        key=lambda x: int(x['count']),
        reverse=True)
    rows = sorted(
        rows,
        key=lambda x: (
            STEP_ORDER.index(x['step']), x['direction'], x['orientation']))
    out.writeheader()
    out.writerows(rows)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
