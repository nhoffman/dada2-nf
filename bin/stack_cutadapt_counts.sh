#!/bin/bash

# Concatenate (stack) cutadapt .tsv summary files,  add sampleid as column

# Includes the header line from the first file; others are
# omitted. Sampleid is extracted from filename which is expected
# as "{sampleid}.cutadapt.tsv"

set -e

if [[ -z $1 ]]; then
    echo "usage: stack_cutadapt_counts.sh sample1.cutadapt.tsv, [sample2.cutadapt.tsv, ...]"
fi

args=( "$@" )

# extract header from first file, add 'sampleid' column in first position
first="${args[0]}"
cat "$first" | head -n 1 | sed -e 's/^/sampleid\t/' | tr -d '\015'

# parse rows from all files and append 'sampleid' value as first column
# needs 'readlink' to handle nextflow symlink renaming
for f in "${args[@]}"; do
    sampleid=$(basename -s ".cutadapt.tsv" "$(readlink -f "$f")")
    cat "$f" | tail -n+2 | sed -e "s/^/$sampleid\t/" | tr -d '\015'
done
