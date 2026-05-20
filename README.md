<p align="center">
  <img src="assets/dada2-nf-logo.svg" alt="dada2-nf logo" width="520">
</p>

[![Run pipeline tests and deploy to registry](https://github.com/nhoffman/dada2-nf/actions/workflows/test-publish.yml/badge.svg)](https://github.com/nhoffman/dada2-nf/actions/workflows/test-publish.yml)
![Tests passing](https://img.shields.io/badge/tests-passing-brightgreen)

# dada2-nf

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
- supports 16S and ITS configurations
- supports dual-index, single-index, and no-index inputs
- can process both forward and reverse orientations when bidirectional output is
  enabled
- writes sequence variant FASTA files, specimen maps, SV tables, weights, and
  count/yield summaries
- supports Docker, Singularity, and local execution

## Configuration

Pipeline inputs are specified via a params file. Example using
`test/minimal/params.json`:

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

Key parameters:

| Parameter | Description |
|-----------|-------------|
| `sample_information` | CSV file mapping samples to index sequences |
| `fastq_list` | Text file listing FASTQ input paths |
| `index_file_type` | One of `dual`, `single`, or `none` |
| `dada_params` | `data/dada_params_250.json` (250 bp reads) or `data/dada_params_300.json` (300 bp reads) |
| `alignment.strategy` | `cmsearch` for 16S, `vsearch` for ITS or custom libraries |
| `nproc` | Number of CPUs (default: 8) |

See `test/` for additional input file examples and `template.json` for
the full parameter schema.

The version of dada2 used in this project can be determined with:

```
docker run --rm ghcr.io/nhoffman/dada2-nf:latest R -q -e 'packageVersion("dada2")'
```

## Execution

Install the Nextflow binary in this directory:

```
wget -qO- https://get.nextflow.io | bash
```

Three profiles are available. All profiles resume from previous runs by
default (`resume = true`) and limit the local executor to 8 concurrent
tasks.

**`standard`** (default) — runs with Singularity, caches images in
`singularity/`:

```
./nextflow run main.nf -params-file test/minimal/params.json
```

**`docker`** — runs with Docker:

```
./nextflow run main.nf -params-file test/minimal/params.json -profile docker
```

**`local`** — runs without a container engine:

```
./nextflow run main.nf -params-file test/minimal/params.json -profile local
```

Common overrides:

```
./nextflow run main.nf \
  -params-file test/minimal/params.json \
  -profile docker \
  --container dada2-nf:local \
  --work_dir work-test \
  --nproc 4
```

## Testing

Run all tests:

```
./test/run_tests.sh
```

Run selected tests by name:

```
./test/run_tests.sh minimal single/vsearch
```

Available tests:

| Test | Description |
|------|-------------|
| `minimal` | Dual-index, cmsearch 16S filtering |
| `single/cmsearch` | Single-index, cmsearch 16S filtering |
| `single/vsearch` | Single-index, vsearch filtering |
| `single/ngs16s` | Single-index, NGS16S filtering |
| `noindex` | No-index input |
| `its` | ITS configuration |

Each test runs the pipeline against a `test/**/params.json` file and
verifies top-level output files against `test/**/base-files.sha256`.

Pass additional Nextflow arguments after `--`. This is how GitHub Actions
runs tests against the freshly built Docker image:

```
./test/run_tests.sh -- -profile docker --container gha_image
```

## Infernal 16S filtering

The covariance model used for Infernal sequence filtering is obtained
from the Rfam database:

- Model: https://rfam.xfam.org/family/RF00177
- Citation: https://rfam.xfam.org/

## Docker image

The Docker image is hosted in the GitHub container registry:
https://github.com/nhoffman/dada2-nf/pkgs/container/dada2-nf

Singularity can transparently ingest the Docker image and create a
locally-cached copy.
