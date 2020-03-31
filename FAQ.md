# Frequently asked question, tips, and tricks

## AWS Batch

AWS Batch is a system for executing Docker Containers inside AWS ECS instances.
A simple AWS Batch system consists of a Compute Environment, a Job Queue,
a Job Definition and a Job.  Simply, a Compute Environment manages EC2/SPOT
instances, the Job Definition is a registered Docker Container, a Job request
is some defined work to be done using that Job Definition (Docker Container) and
a Job Queue is the gateway between Job requests and the Compute Environment.

https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html

### Simple example

Using this utility:

https://github.com/crosenth/aws_batch

A simple Batch request (with some verbosity) looks like this:

```
aws_batch -v --job-queue optimal --command 'echo hello world; echo and goodbye!' ubuntu-18-04
mkdir -p tmp; cd tmp; echo hello world
ubuntu-18-04 "echo hello world; echo and goodbye!"
RUNNABLE
STARTING
hello world
and goodbye!
SUCCEEDED
```

### What Compute Environments, Queues and Job Definitions are available?

Using the awscli (pip install awscli):

```
aws batch describe-compute-environments
...
aws batch describe-queues
...
aws batch describe-job-definitions
...
```

For a full list of commands:

```
aws batch help
```

Note: many of the commands are only available to AWS admins

### How do I understand resource allocation?

Resource allocation defines which EC2/SPOT instance will be chosen to run your
Job request. Resource allocation can be defined at any point in the Batch
process with the Batch Job having the final word on how much compute resources
it needs to execute properly.  EC2 instances come in a variety of CPU and RAM
combinations:

https://aws.amazon.com/ec2/instance-types/

The minimum requirements for any Batch Job are 1 cpu and 
4 MiB of RAM.

### How do I create my own Job Definitions?

Docker containers are registered with the following command:

```
aws batch register-job-definition --cli-input-json file://jobDefinition.json
```

The jobDefinition.json file could look something like this:

```
% cat jobDefininition.json
{
    "jobDefinitionName": "myjobDefinition",
    "type": "container",
    "containerProperties": {
        "image": "ubuntu:18.04",
        "vcpus": 1,
        "memory": 1024,
        "command": [ "true" ],
        "volumes": [
            {
                "host": {
                    "sourcePath": "/home/ec2-user/miniconda"
                },
                "name": "aws-cli"
            },
        ],
				"mountPoints": [
            {
                "containerPath": "/home/ec2-user/miniconda",
                "readOnly": true,
                "sourceVolume": "aws-cli"
            },
        ],
    }
}
```

In this example the awscli is mounted from the EC2 instance into the Docker
Container to the give users the ability to move data in and out of the
Container for processing.  The awscli can also be installed to the Container
itself from Docker Image and Dockerfile.

Also note: Containers without an entry point must be configured
`"command": [ "true"  ].

### Specifying job resources

If a Job is expected to exceed the resource requirements defined in a Job
Definition the Job can override the Job Definition.
Defining a higher number of CPUs and more memory will also be
mapped to the Docker run arguments `--cpu-shares` and `--memory`:

https://docs.docker.com/config/containers/resource_constraints/

Keep in mind `--cpu-shares` will not prevent Docker containers from utilizing
every CPU available in an EC2 instance but `--memory` *will* prevent a Docker
container from utilizing more than the requested amount.  Think of CPU
allocation as a floor requirment and the RAM requirement as a ceiling
requirement.

### How can I tell which EC2 instance type was allocated for my Job?

The EC2 instance can be determined using the AWS Instance Metadata Service by
executing the following using wget (or curl) as part of the Job request:

```
aws_batch --job-queue optimal --command "wget -O - -q http://169.254.169.254/latest/meta-data/instance-type" ubuntu-18-04
c4.2xlarge
```

For more information on using the Instance Metadata Service see

https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html

### Can I view other Job requests being executed on a container instance?

You can view all Job tasks currently being executed within an EC2 Instance
using AWS ECS Container Agent Introspection by running wget (or curl) on
port 51678 from localhost:

```
wget -O - -q http://localhost:51678/v1/tasks"
```

For more information on using ECS Container Agent Introspection:

https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-agent-introspection.html

### How can I test my aws credentials?

Easiest way to test or generate aws credentials is to execute the following
using the awscli tool (pip install awscli):

```
aws configure
```

By default the credentials are stored in the $HOME directory in the
`~/.aws/credentials` file.  A further test of credentials could be attempting
to view an AWS S3 bucket (see below).

### How can I see files in my s3 bucket folder?

```
 aws s3 ls --human-readable s3://mybucket/folder/
```

### How do I clean up my s3 bucket folder?

Use the `--recursive` command to remove an entire folder:

```
aws s3 rm --recursive s3://mybucket/folder/
```

Or if you just want to remove a single file:

```
aws s3 rm s3://mybucket/folder/remove-me.txt
```

For information on additional, useful options:

```
aws s3 rm help
```

## Pipeline development

### How do I force local execution of a task without Docker?

Include ``container null`` within the process definition.

### How can I test and develop individual steps of the pipeline?

This pipeline was designed to minimize the logic and inline code in
Nextflow. Invocations of the ``dada2`` package are make via scripts in
``bin/``. One pattern for testing and development is to execute an
individual script in the docker environment using files produced by
local execution of the pipeline. For example:

```
docker run --rm -it -v $(pwd):$(pwd) -w $(pwd) \
	quay.io/nhoffman/dada2-nf:v1.12-9-g39fa45b bin/dada2_learn_errors.R \
	output/batch_1/filtered --model errors.rds
```

This makes it possible to iteratively develop the scripts, add print
statements, examine output, etc, outside the context of the pipeline.

### What does this Nextflow error mean?

```
Process `...` terminated for an unknown reason -- Likely it has been terminated by the external system
```

This error indicates a problem with the Job Definition configuration or the
location of the awscli tool.

Also see: https://www.nextflow.io/docs/latest/awscloud.html#troubleshooting


