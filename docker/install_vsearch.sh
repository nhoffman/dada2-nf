#!/bin/bash

# Install vsearch binaries to $prefix/bin

# https://github.com/torognes/vsearch/releases/download/v2.17.0/vsearch-2.17.0-linux-aarch64.tar.gz

VERSION=2.17.0
VSEARCH=vsearch-${VERSION}-linux-x86_64

cd /tmp
wget https://github.com/torognes/vsearch/releases/download/v${VERSION}/${VSEARCH}.tar.gz
tar xzf "${VSEARCH}.tar.gz"
cp ${VSEARCH}/bin/* "/usr/local/bin"
rm -rf "${VSEARCH}"
