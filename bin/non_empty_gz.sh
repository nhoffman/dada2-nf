#!/bin/bash

for file in "$@"
do
    # https://unix.stackexchange.com/a/6771
    if LC_ALL=C gzip -l "$file" | awk 'NR==2 {exit($2!=0)}'; then
	:
    else
	echo "$file"
    fi
done
