#!/usr/bin/env python3
"""Check automatic manifest paths and all supported input branches."""

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class DirectoryInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='dada-dir-test-')
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.reads = self.path / 'reads with spaces'
        self.reads.mkdir()
        for i, sample in enumerate(('alpha', 'beta'), 1):
            for mate in (1, 2):
                name = f'{sample}_S{i}_R{mate}_001.fastq.gz'
                with gzip.open(self.reads / name, 'wt') as handle:
                    handle.write(f'@{sample}\nACGT\n+\nIIII\n')

    def make_manifest(self, source=None, cwd=None):
        output = self.path / 'manifest.csv'
        subprocess.run([
            sys.executable, str(ROOT / 'bin' / 'create_manifest.py'),
            str(source or self.reads), '--batch', '1', '--out', str(output)
        ], cwd=cwd or self.path, check=True, capture_output=True)
        return output

    def preview(self, inputs):
        if not shutil.which('nextflow'):
            self.skipTest('Nextflow is required for workflow previews')
        params = {
            'dada_params': str(ROOT / 'data' / 'dada_params_300.json'),
            'index_file_type': 'none',
            'alignment': {'strategy': 'none'},
            'output': str(self.path / 'output'),
            'work_dir': str(self.path / 'work'),
            'nproc': 1,
            'bidirectional': False,
        }
        params.update(inputs)
        config = self.path / 'params.json'
        config.write_text(json.dumps(params))
        env = dict(os.environ, NXF_DISABLE_CHECK_LATEST='true')
        result = subprocess.run([
            'nextflow', '-log', str(self.path / 'nextflow.log'), 'run',
            str(ROOT / 'main.nf'), '-params-file', str(config),
            '-profile', 'local', '-preview', '-w', str(self.path / 'work'),
            '-ansi-log', 'false', '-offline'
        ], cwd=self.path, env=env, capture_output=True, text=True,
           timeout=60)
        return result, result.stdout + result.stderr

    def make_sample_sheet(self):
        sheet = self.path / 'sample_sheet.csv'
        with sheet.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sampleid', 'batch'])
            writer.writerows([['alpha', 7], ['beta', 8]])
        return sheet

    def make_fastq_list(self):
        paths = self.path / 'fastq_list.txt'
        paths.write_text(''.join(str(path) + '\n' for path in
                                 sorted(self.reads.glob('*.fastq.gz'))))
        return paths

    def test_staged_directory_keeps_original_paths_and_batch_one(self):
        work = self.path / 'task'
        work.mkdir()
        staged = work / 'fastq_input'
        staged.symlink_to(self.reads, target_is_directory=True)
        output = self.make_manifest(staged, work)
        with output.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row['sampleid'] for row in rows},
                         {'alpha', 'beta'})
        for row in rows:
            self.assertEqual(row['batch'], '1')
            self.assertEqual(row['datadir'], str(self.reads.resolve()))
            for mate in ('R1', 'R2'):
                self.assertTrue((Path(row['datadir']) / row[mate]).is_file())
            self.assertEqual(row['I1'], '')
            self.assertEqual(row['I2'], '')

    def test_directory_is_the_only_input_selector(self):
        result, text = self.preview({'fastq_dir': str(self.reads)})
        self.assertEqual(result.returncode, 0, text)

    def test_manifest_takes_precedence_over_other_inputs(self):
        result, text = self.preview({
            'manifest': str(self.make_manifest()),
            'fastq_dir': str(self.path / 'missing-directory'),
            'sample_information': str(self.path / 'missing-samples.csv'),
            'sample_sheet': str(self.path / 'missing-sheet.csv'),
            'fastq_list': str(self.path / 'missing-fastqs.txt'),
        })
        self.assertEqual(result.returncode, 0, text)

    def test_empty_directory_is_rejected(self):
        empty = self.path / 'empty'
        empty.mkdir()
        result, text = self.preview({'fastq_dir': str(empty)})
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn('No .fastq.gz files found', text)

    def test_regular_file_is_rejected(self):
        read = next(self.reads.glob('*.fastq.gz'))
        result, text = self.preview({'fastq_dir': str(read)})
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn('--fastq_dir must point to a directory', text)

    def test_missing_directory_is_rejected(self):
        missing = self.path / 'missing'
        result, text = self.preview({'fastq_dir': str(missing)})
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn(str(missing), text)

    def test_existing_manifest_mode(self):
        result, text = self.preview({'manifest': str(self.make_manifest())})
        self.assertEqual(result.returncode, 0, text)

    def test_existing_sample_information_and_fastq_list_mode(self):
        manifest = self.make_sample_sheet()
        paths = self.make_fastq_list()
        result, text = self.preview({
            'sample_information': str(manifest), 'fastq_list': str(paths)
        })
        self.assertEqual(result.returncode, 0, text)

    def test_sample_information_and_list_take_precedence_over_directory(self):
        manifest = self.make_sample_sheet()
        paths = self.make_fastq_list()
        result, text = self.preview({
            'sample_information': str(manifest), 'fastq_list': str(paths),
            'fastq_dir': str(self.path / 'missing-directory'),
        })
        self.assertEqual(result.returncode, 0, text)

    def test_sample_sheet_with_fastq_list(self):
        result, text = self.preview({
            'sample_sheet': str(self.make_sample_sheet()),
            'fastq_list': str(self.make_fastq_list()),
        })
        self.assertEqual(result.returncode, 0, text)

    def test_sample_sheet_with_directory(self):
        sheet = self.make_sample_sheet()
        before = sheet.read_bytes()
        result, text = self.preview({
            'sample_sheet': str(sheet), 'fastq_dir': str(self.reads),
        })
        self.assertEqual(result.returncode, 0, text)
        self.assertEqual(sheet.read_bytes(), before)

    def test_legacy_sample_information_with_directory(self):
        result, text = self.preview({
            'sample_information': str(self.make_sample_sheet()),
            'fastq_dir': str(self.reads),
        })
        self.assertEqual(result.returncode, 0, text)

    def test_sample_sheet_requires_a_fastq_source(self):
        result, text = self.preview({
            'sample_sheet': str(self.make_sample_sheet()),
        })
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn('A sample sheet requires fastq_list or fastq_dir', text)

    def test_explicit_list_takes_precedence_with_a_sample_sheet(self):
        result, text = self.preview({
            'sample_sheet': str(self.make_sample_sheet()),
            'fastq_list': str(self.make_fastq_list()),
            'fastq_dir': str(self.path / 'missing-directory'),
        })
        self.assertEqual(result.returncode, 0, text)

    def test_no_input_is_rejected(self):
        result, text = self.preview({})
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn('Provide manifest, sample_sheet/sample_information',
                      text)


if __name__ == '__main__':
    unittest.main()
