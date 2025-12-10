FROM rocker/r-ubuntu:24.04

ARG DADA2_REF=v1.26
ENV DADA2_REF=$DADA2_REF PIP_NO_CACHE_DIR=1

# https://ethicalhackx.com/speed-apt-get-update-parallel-downloads/
ADD 99parallel /etc/apt/apt.conf.d/99parallel
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install --no-install-recommends -y \
    fasttree \
    python3-dev\
    r-cran-ape \
    r-cran-biocmanager \
    r-cran-devtools \
    r-cran-dplyr \
    r-cran-lattice \
    r-cran-latticeextra \
    r-cran-gridextra \
    r-cran-r.utils \
    r-cran-readr \
    r-cran-rmarkdown \
    r-cran-tidyr \
    libcurl4-openssl-dev \
    libxml2 \
    libxml2-dev \ 
    wget && \
    rm -rf /var/lib/apt/lists/*

ADD install-dada2.R /tmp/
# uses value of $DADA2_REF
RUN Rscript /tmp/install-dada2.R
RUN mkdir -p /app /fh /mnt /run/shm

ADD install_infernal_and_easel.sh /tmp/
RUN /tmp/install_infernal_and_easel.sh

ADD install_vsearch.sh /tmp/
RUN /tmp/install_vsearch.sh 2.30.1

ADD install_xsv.sh /tmp/
RUN /tmp/install_xsv.sh

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

ADD requirements.txt /tmp/
ENV VIRTUAL_ENV=/usr/local/venv
RUN uv venv $VIRTUAL_ENV
RUN uv pip install -r /tmp/requirements.txt

ENV PATH="$VIRTUAL_ENV/bin:$PATH"
