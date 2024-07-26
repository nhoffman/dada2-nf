FROM rocker/r-ubuntu:22.04

ARG DADA2_REF=v1.18
ENV DADA2_REF=$DADA2_REF

# https://ethicalhackx.com/speed-apt-get-update-parallel-downloads/
ADD 99parallel /etc/apt/apt.conf.d/99parallel
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
    fasttree \
    python3-pip
#     r-cran-ape \
#     r-cran-biocmanager \
#     r-cran-devtools \
#     r-cran-dplyr \
#     r-cran-lattice \
#     r-cran-latticeextra \
#     r-cran-gridextra \
#     r-cran-r.utils \
#     r-cran-readr \
#     r-cran-rmarkdown \
#     r-cran-tidyr \
#     libcurl4-openssl-dev \
#     libxml2 \
#     libxml2-dev \
#     vsearch
# 
# ADD install-dada2.R /tmp/
# # uses value of $DADA2_REF
# RUN Rscript /tmp/install-dada2.R
# RUN mkdir -p /app /fh /mnt /run/shm
# 
# ADD install_infernal_and_easel.sh /tmp/
# RUN /tmp/install_infernal_and_easel.sh
# 
# ADD install_vsearch.sh /tmp/
# RUN /tmp/install_vsearch.sh 2.21.1
# 
# ADD install_xsv.sh /tmp/
# RUN /tmp/install_xsv.sh
# 
# ADD install_pplacer.sh /tmp/
# RUN /tmp/install_pplacer.sh

ADD requirements.txt /tmp/
RUN pip3 install -U pip wheel
RUN pip3 install -r /tmp/requirements.txt

ADD test/ /tmp/
