FROM rocker/r-apt:bionic

ARG DADA2_COMMIT
ENV DADA2_COMMIT=$DADA2_COMMIT

ADD 99parallel /etc/apt/apt.conf.d/99parallel
RUN apt-get update -y && apt-get install -y \
    python3-pip \
    r-cran-ape \
    r-cran-biocmanager \
    r-cran-devtools \
    r-cran-dplyr \
    r-cran-lattice \
    r-cran-latticeextra \
    r-cran-r.utils \
    r-cran-readr \
    r-cran-rmarkdown \
    r-cran-tidyr \
    libcurl4-openssl-dev \
    libxml2 \
    libxml2-dev

ADD requirements.txt /tmp/
RUN pip3 install -r /tmp/requirements.txt

ADD install_infernal_and_easel.sh /tmp/
RUN /tmp/install_infernal_and_easel.sh

ADD install-dada2.R /tmp/
RUN Rscript /tmp/install-dada2.R
RUN mkdir -p /app /fh /mnt /run/shm

