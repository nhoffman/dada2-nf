<p align="center">
  <img src="assets/dada2-nf-logo.svg" alt="dada2-nf logo" width="520">
</p>

[![Run pipeline tests and deploy to registry](https://github.com/nhoffman/dada2-nf/actions/workflows/test-publish.yml/badge.svg)](https://github.com/nhoffman/dada2-nf/actions/workflows/test-publish.yml)
![Tests passing](https://img.shields.io/badge/tests-passing-brightgreen)

# Dada2 Nextflow pipeline

Filters, trims, and denoises NGS short reads using
[dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html).
Additional features:

- optionally removes adapters with [cutadapt](https://cutadapt.readthedocs.io)
- filters indexed reads with [barcodecop](https://github.com/nhoffman/barcodecop)
- generates quality plots for forward and reverse reads
- learns DADA2 error models, denoises reads, merges paired reads, and removes
  chimeras
- filters off-target reads using cmsearch with covariance models or vsearch
  with reference libraries
- supports 16S and ITS test configurations
- supports dual-index, single-index, and no-index inputs
- can process both forward and reverse orientations when bidirectional output is
  enabled
- writes sequence variant FASTA files, specimen maps, SV tables, weights, and
  count/yield summaries
- supports Docker, Singularity, local execution, S3 inputs, and AWS Batch
  configuration

Example configuration:

```json
{
  "sample_information": "test/minimal/sample-information.csv",
  "fastq_list": "test/minimal/fastq-list.txt",
  "manifest": "",
  "output": "output-minimal",
  "index_file_type": "dual",
  "dada_params": "data/dada_params_300.json",
  "bidirectional": false,
  "alignment": {
    "library": "",
    "model": "data/SSU_rRNA_bacteria.cm",
    "strategy": "cmsearch"
  }
}
```

See contents of the ``test/`` directory for examples of input files.

The version of dada2 used in this project (which has no relation to
the version tag for this repository) can be determined with this command:

```
% docker run --rm -it ghcr.io/nhoffman/dada2-nf:latest R -q -e 'packageVersion("dada2")'
> packageVersion("dada2")
[1] ‘1.18.0’
```

## Execution options

Install the Nextflow binary in this directory:

```
wget -qO- https://get.nextflow.io | bash
```

The default profile is ``standard``. It runs local processes using Singularity
using the configured container image, caches Singularity images in
``singularity/`` and limits the local executor queue to 8 concurrent tasks.

```
./nextflow run main.nf -params-file test/minimal/params.json
```

To run with Docker instead of Singularity, use the ``docker`` profile:

```
./nextflow run main.nf -params-file test/minimal/params.json -profile docker
```

To run without enabling a container engine, use the ``local`` profile:

```
./nextflow run main.nf -params-file test/minimal/params.json -profile local
```

Runtime settings from ``nextflow.config`` can be overridden on the command
line. Common overrides include:

```
./nextflow run main.nf \
  -params-file test/minimal/params.json \
  -profile docker \
  --container dada2-nf:local \
  --work_dir work-test \
  --nproc 4
```

## Testing

Run all regression tests:

```
./test/run_tests.sh
```

Run selected tests by name:

```
./test/run_tests.sh minimal single/vsearch
```

Pass additional Nextflow arguments after ``--``. This is how the GitHub Actions
workflow runs the checksum tests against the freshly built Docker image:

```
./test/run_tests.sh -- -profile docker --container gha_image
```

Each test runs the pipeline with a ``test/**/params.json`` file and then checks
top-level deterministic output files against the matching
``test/**/base-files.sha256`` file.

## Infernal 16S filtering

Covariance model used for Infernal sequence filtering obtained from the Rfam database:

https://rfam.xfam.org/family/RF00177

To cite Rfam see latest web site instructions:

https://rfam.xfam.org/

## Docker image

A Docker image is hosted in the GitHub container registry:
https://github.com/nhoffman/dada2-nf/pkgs/container/dada2-nf

Singularity can transparently ingest the Docker image and create a locally-cached image.
