#!/bin/bash

set -e

# if ! git diff-index --quiet HEAD; then
#     echo "Error: git repo is dirty - commit and try again"
#     echo
#     git status
#     exit 1
# fi

repo=dada2-nf
dada2_version=v${1-1.18}
rev=${2-$(git describe --tags --dirty)}

python3 get_tag.py $dada2_version > /dev/null
DADA2_COMMIT=$(python3 get_tag.py $dada2_version)
image="${repo}:v${rev}"

echo "building image $image from dada2 $dada2_version commit $DADA2_COMMIT"
docker build --build-arg DADA2_COMMIT=$DADA2_COMMIT --rm --force-rm -t "$image" .

