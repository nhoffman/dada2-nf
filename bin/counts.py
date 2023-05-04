#!/usr/bin/env python3
import argparse
import csv
import sys

STEP_ORDER = [
    'raw',
    'barcodecop',
    'downsample',
    'cutadapt',
    'on_target',
    'filter_and_trim',
    'dada2_denoise',
    'dada2_merge',
    'dada2_chimera',
    'svs']

DIR_ORDER = ['R1', 'R2', 'merged', '']


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest', type=argparse.FileType('r'))
    parser.add_argument('downsample', type=int)
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

    # raw
    raw = list(csv.DictReader(args.manifest))
    yld = {r['sampleid']: int(r['count']) for r in raw}
    rows.extend(process_rows(raw, yld, 'raw', 'count'))

    # barcodecop
    barcodecop = csv.DictReader(
        args.barcodecop,
        fieldnames=['filename', 'in', 'count'])
    bcop = {}
    # deduplicate barcodecop rows by sampleid from filename
    for b in barcodecop:
        si = b['filename'].split('_')[0]
        b['sampleid'] = si
        bcop[si] = b
    bcop = bcop.values()
    bcop = list(process_rows(bcop, yld, 'barcodecop', 'count'))
    rows.extend(bcop)

    if args.downsample != -1:
        for row in bcop:
            si = row['sampleid']
            if args.downsample < yld[si]:
                yld[si] = args.downsample

    # cutadapt
    cutadapt = list(csv.DictReader(args.cutadapt))
    if args.downsample != -1:
        rows.extend(process_rows(cutadapt, yld, 'downsample', 'in_reads'))
    rows.extend(process_rows(cutadapt, yld, 'cutadapt', 'out_reads'))

    # on_target
    splits = list(csv.DictReader(
        args.split_orientations,
        fieldnames=['sampleid', 'orientation', 'count']))
    rows.extend(process_rows(splits, yld, 'on_target', 'count'))

    # dada2
    dada2 = list(csv.DictReader(args.dada2))
    rows.extend(
        process_rows(dada2, yld, 'filter_and_trim', 'filtered_and_trimmed'))
    rows.extend(
        process_rows(dada2, yld, 'dada2_denoise', 'denoised_r1', 'R1'))
    rows.extend(
        process_rows(dada2, yld, 'dada2_denoise', 'denoised_r2', 'R2'))
    rows.extend(
        process_rows(dada2, yld, 'dada2_merge', 'no_chimeras', 'merged'))
    chimeras = [d for d in dada2 if d['merged'] != d['no_chimeras']]
    for c in chimeras:
        c['chimera'] = int(c['merged']) - int(c['no_chimeras'])
    rows.extend(
        process_rows(chimeras, yld, 'dada2_chimera', 'chimera', 'merged'))

    # svs
    fieldnames = ['sampleid', 'direction', 'count']
    svs = list(csv.DictReader(args.specimens, fieldnames=fieldnames))
    merged = [s for s in svs if s['direction'] == 'merged']
    rows.extend(process_rows(merged, yld, 'svs', 'count', 'merged'))
    r1 = [s for s in svs if s['direction'] == 'R1']
    rows.extend(process_rows(r1, yld, 'svs', 'count', 'R1'))
    r2 = [s for s in svs if s['direction'] == 'R2']
    rows.extend(process_rows(r2, yld, 'svs', 'count', 'R2'))

    # output
    out = csv.DictWriter(
        args.out,
        fieldnames=[
            'step',
            'sampleid',
            'orientation',
            'direction',
            'count',
            'denominator',
            'yield'],
        extrasaction='ignore')
    rows = sorted(
        rows,
        key=lambda x: int(x['count']),
        reverse=True)
    rows = sorted(
        rows,
        key=lambda x: (
            STEP_ORDER.index(x['step']),
            DIR_ORDER.index(x['direction']),
            x['orientation']))
    out.writeheader()
    out.writerows(rows)


def process_rows(rows, yld, step, count, direction=''):
    for r in rows:
        si = r['sampleid']
        yield {
            'count': r[count],
            'direction': direction,
            'orientation': '',
            'step': step,
            'yield': int(r[count]) / yld[si] if yld[si] else 0,
            'denominator': yld[si],
            **r}


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
