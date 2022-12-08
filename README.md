# Dada2 Nextflow pipeline

Filters, trims, and denoises NGS short reads using
[dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html). Additional features:

- optionally removes adapters with [cutadapt](https://cutadapt.readthedocs.io)
- filters reads with [barcodecop](https://github.com/nhoffman/barcodecop)
- generates quality plots for forward and reverse reads
- provides a table of sequence yields at each step of the pipeline
- removes non-16S rRNA gene sequences uses cmsearch and an Rfam alignment model
- compatibility with S3 objects as input and AWS Batch

Example configuration:

```
% cat params-minimal.json
{
  "sample_information": "test/sample-information-minimal.csv",
  "fastq_list": "test/fastq-list-minimal.txt",
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

## Local execution quickstart for the truly impatient

Install the nextflow binary in this directory

```
wget -qO- https://get.nextflow.io | bash
```

Execute locally with the default Docker image using the minimal data set.

```
./nextflow run main.nf -params-file params-minimal.json
```

Local execution using the Singularity image defined in ``nextflow.config``

```
./nextflow run main.nf -params-file params-minimal.json -profile singularity
```

An alternative Docker or Singularity image (eg, a version other than
``:latest`` or one that is local) may be specified in ``-params-file``
by adding a "container" element, or as an argument to the command line
argument ``--container``.

Profiles that run locally (see ``nextflow.config``) use a default
workDir named "work"; another name can be specified using the command
line argument ``--work_dir``.

## Execution on AWS Batch

Details will depend on your AWS batch configuration. General instructions TBD.

## Infernal 16s filtering

Coveriance model used for Infernal sequence filtering obtained from the Rfam database:

https://rfam.xfam.org/family/RF00177

To cite Rfam see latest web site instructions:

https://rfam.xfam.org/

## Docker image

A Docker image is hosted in the GitHub container registry:
https://github.com/nhoffman/dada2-nf/pkgs/container/dada2-nf

Singularity can transparently ingest the Docker image and create a locally-cached image.
