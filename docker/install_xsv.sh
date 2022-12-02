#!/bin/bash

# Install xsv binaries to $prefix/bin
VERSION=0.13.0
XSV=xsv-$VERSION-x86_64-unknown-linux-musl
URL=https://github.com/BurntSushi/xsv/releases/download/$VERSION/$XSV.tar.gz
echo $URL
cd /tmp
wget --quiet -nc $URL
tar -xf "${XSV}.tar.gz"
cp xsv "/usr/local/bin"
rm -rf "${XSV}"
