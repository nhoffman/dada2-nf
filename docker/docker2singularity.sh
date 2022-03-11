#!/bin/bash

last_tag=$(docker images dada2-nf --format '{{.Tag}}' | head -n1)
tag=${1-$last_tag}
echo "Buiding Singularity image for dada2-nf:$tag"

docker run \
       -v "${DOCKER_HOST#unix://}:/var/run/docker.sock" \
       -v "$(pwd):/output" --privileged -t --rm \
       quay.io/singularity/docker2singularity \
       --name dada2-nf-${tag}.sif dada2-nf:$tag

