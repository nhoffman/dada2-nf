#!/bin/bash

# iterate over the filenames of provided fastq.gz
# and print the filenames of only non-empty files
for file in "$@"
do
    # check the gzip -l output 'Uncompressed' size
    # and if zero, exit. If non-zero, return the filename to stdout
    # See https://unix.stackexchange.com/a/6771 for the hint
    if LC_ALL=C gzip -l "$file" | awk 'NR==2 {exit($2==0)}'; then
	echo "$file"
    fi
done
