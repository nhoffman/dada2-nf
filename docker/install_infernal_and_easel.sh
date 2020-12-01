#!/bin/bash

# Install infernal (cmalign) and selected easel binaries to $prefix/bin

# http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz

VERSION=1.1.2
INFERNAL=infernal-${VERSION}-linux-intel-gcc
cd /tmp
wget -q -nc http://eddylab.org/infernal/${INFERNAL}.tar.gz
for binary in cmalign cmconvert cmsearch esl-alimerge esl-sfetch esl-reformat; do
    tar xf "${INFERNAL}.tar.gz" --no-anchored "binaries/$binary"
done
cp ${INFERNAL}/binaries/* "/usr/local/bin"
rm -r "${INFERNAL}"


