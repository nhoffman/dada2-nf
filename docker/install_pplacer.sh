#!/bin/bash

bin=/usr/local/bin

# install pplacer and accompanying python scripts
PPLACER_BUILD=1.1.alpha19
PPLACER_DIR=pplacer-Linux-v${PPLACER_BUILD}
PPLACER_ZIP=${PPLACER_DIR}.zip
PPLACER_GH="https://github.com/matsen/pplacer"

pplacer_is_installed(){
    $bin/pplacer --version 2> /dev/null | grep -q "$PPLACER_BUILD"
}

if pplacer_is_installed; then
    echo -n "pplacer is already installed: "
    $bin/pplacer --version
else
    mkdir -p src && \
	(cd src && \
	wget -nc --quiet $PPLACER_GH/releases/download/v$PPLACER_BUILD/$PPLACER_ZIP && \
	unzip -o $PPLACER_ZIP && \
	cp $PPLACER_DIR/{pplacer,guppy,rppr} $bin)
    # confirm that we have installed the requested build
    if ! pplacer_is_installed; then
	echo -n "Error: you requested pplacer build $PPLACER_BUILD "
	echo "but $($venv/bin/pplacer --version) was installed."
	echo "Try removing src/$PPLACER_ZIP first."
	exit 1
    fi
fi

