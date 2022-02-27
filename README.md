# Dada2 Nextflow pipeline

Filters, trims, and denoises NGS short reads using [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html). Additional features:

- filters reads with [barcodecop](https://github.com/nhoffman/barcodecop)
- generates quality plots for forward and reverse reads
- provides a table of sequence yields at each step of the pipeline
- removes non-16S rRNA gene sequences uses cmsearch and an Rfam alignment model
- compatibility with S3 objects as input

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
  "alignment_model": "data/SSU_rRNA_bacteria.cm"
}
```

See contents of the ``test/`` directory for examples of input files.

## Local execution quickstart for the truly impatient

Install the nextflow binary in this directory

```
wget -qO- https://get.nextflow.io | bash
```

Execute locally with Docker, using the minimal data set.

```
./nextflow run main.nf -params-file params-minimal.json
```

Local execution using the Singularity image defined in ``nextflow.config``

```
./nextflow run main.nf -params-file params-minimal.json -profile singularity
```

An alternative Singularity image (eg, one that is local) may be specified in ``-params-file`` , eg:

```
  "singularity_container": "dada2-nf_v1.15-dev-2022-02-15-a68e40f4dd5b.sif"
```

## Execution on AWS Batch

Details will depend on your AWS batch configuration. General instructions TBD.

### Infernal 16s filtering

Coveriance model used for Infernal sequence filtering obtained from the Rfam database:

https://rfam.xfam.org/family/RF00177

To cite Rfam see latest web site instructions:

https://rfam.xfam.org/

## Docker and Singularity images

- Docker Hub: https://hub.docker.com/repository/docker/nghoffman/dada2-nf
- Singularity (Sylabs.io): https://cloud.sylabs.io/library/nhoffman/dada2-nf/dada2-nf
