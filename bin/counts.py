#!/usr/bin/env python3
import argparse
import csv
import sys

ORIENTATIONS_ORDER = ['forward', 'reverse', 'off_target']


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', type=argparse.FileType('r'))
    parser.add_argument('cutadapt', type=argparse.FileType('r'))
    parser.add_argument('split_orientations', type=argparse.FileType('r'))
    parser.add_argument('barcodecop', type=argparse.FileType('r'))
    parser.add_argument('dada2', type=argparse.FileType('r'))
    parser.add_argument('passed', type=argparse.FileType('r'))
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = parser.parse_args(arguments)
    raw = {r['sampleid']: r['count'] for r in csv.DictReader(args.manifest)}
    cutadapt = csv.DictReader(args.cutadapt)
    cutadapt = (c for c in cutadapt if c['sampleid'] in raw)
    cutadapt = {c['sampleid']: c for c in cutadapt}
    splits = csv.DictReader(
        args.split_orientations,
        fieldnames=['sampleid', 'orientation', 'reoriented'])
    splits = (s for s in splits if s['sampleid'] in raw)
    splits = {(s['sampleid'], s['orientation']): s for s in splits}
    barcodecop = csv.reader(args.barcodecop)
    barcodecop = ((b[0].split('_')[0], b[2]) for b in barcodecop)
    barcodecop = (b for b in barcodecop if b[0] in raw)
    barcodecop = dict(barcodecop)
    dada2 = csv.DictReader(args.dada2)
    dada2 = (d for d in dada2 if d['sampleid'] in raw)
    dada2 = {(d['sampleid'], d['orientation']): d for d in dada2}
    passed = csv.DictReader(args.passed)
    passed = (p for p in passed if p['sampleid'] in raw)
    passed = {(p['sampleid'], p['orientation']): p for p in passed}
    out = csv.DictWriter(
        args.out,
        fieldnames=[
            'sampleid',
            'raw',
            'barcodecop',
            'cutadapt',
            'orientation',
            'reoriented',
            'filtered_and_trimmed',
            'denoised_r1',
            'denoised_r2',
            'merged',
            'no_chimeras',
            'target',
            'not_target'],
        extrasaction='ignore')
    out.writeheader()
    samples = sorted(
        splits.keys() | dada2.keys() | passed.keys(),
        key=lambda x: (x[0], ORIENTATIONS_ORDER.index(x[1])))
    for sampleid, orientation in samples:
        if sampleid in cutadapt and 'out_reads' in cutadapt[sampleid]:
            cuta = cutadapt[sampleid]['out_reads']
        else:
            cuta = ''
        out.writerow({
            'sampleid': sampleid,
            'raw': raw[sampleid],
            'barcodecop': barcodecop[sampleid],
            'cutadapt': cuta,
            'orientation': orientation,
            **splits[(sampleid, orientation)],
            **dada2.get((sampleid, orientation), {}),
            **passed.get((sampleid, orientation), {})})


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
