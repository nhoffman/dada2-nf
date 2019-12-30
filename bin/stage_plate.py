#!/usr/bin/env python3

import sys
import argparse
from os import path
from glob import glob
import json
import os
import subprocess


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('plate', help="plate label, ie miseq-plate-{label}")
    parser.add_argument('-d', '--data-dir',
                        default='/fh/fast/fredricks_d/bvdiversity/data')

    args = parser.parse_args(arguments)
    plate = args.plate

    platedir = path.join(args.data_dir, f'miseq-plate-{plate}')
    outdir = path.join(args.data_dir, 'dada2_nf_out', f'miseq-plate-{plate}')
    try:
        os.makedirs(outdir)
    except OSError as err:
        pass

    # fastq list
    fq_pattern = path.join(
        platedir,
        f'run-files/*/Data/Intensities/BaseCalls/m{plate}*.fastq.gz')
    files = glob(fq_pattern)

    fastq_list = path.abspath(path.join(outdir, 'fastq_list.txt'))
    with open(fastq_list, 'w') as f:
        f.write('\n'.join(files) + '\n')

    # sample information
    sample_information = path.join(
        platedir, f'sample-information/sample-information-m{plate}.xlsx')
    assert os.path.exists(sample_information)

    params_file = path.join(outdir, 'params.json')
    d = {
        'sample_information': sample_information,
        'fastq_list': fastq_list,
        'output': path.join(outdir, 'output'),
        'index_file_type': 'single',
        'dada_params': 'data/dada_params_300.json',
    }
    with open(params_file, 'w') as f:
        json.dump(d, f, indent=2)

    version_info = path.join(outdir, 'version_info.txt')
    cmds = [
        f'date > {version_info}',
        f'pwd >> {version_info}',
        f'git describe --tags --dirty >> {version_info}',
        f'git log -n1 >> {version_info}',
    ]
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

    print(f'nextflow run main.nf -profile hutch_batch -params-file {params_file}')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
