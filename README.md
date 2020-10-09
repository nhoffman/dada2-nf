# Dada2 Nextflow pipeline

## Local execution quickstart for the truly impatient

Install Docker and make sure that the Docker daemon is running.

Install the nextflow binary in this directory

```
wget -qO- https://get.nextflow.io | bash
```

Execute locally, using the minimal data set.

```
./nextflow run main.nf -params-file params-minimal.json
```

## Execution on AWS Batch

Details will depend on your AWS batch configuration. General instructions TBD.

### Infernal 16s filtering

Coveriance model used for Infernal sequence filtering obtained from the Rfam database:

https://rfam.xfam.org/family/RF00177

To cite Rfam see latest web site instructions:

https://rfam.xfam.org/
