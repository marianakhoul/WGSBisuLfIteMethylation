FROM ubuntu:20.04

RUN apt update \
    && apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common -y \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    #&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt install pandoc curl libcurl4-openssl-dev libssl-dev libxml2-dev r-base r-base-core r-recommended -y #\
    #&& Rscript -e "install.packages(c('Rserve', 'data.table', 'BiocManager', 'KernSmooth', 'parallel', 'tidyverse', 'shiny'))" \
    #&& Rscript -e "BiocManager::install(c('GenomicRanges', 'BSgenome', 'bsseq', 'MethylSeekR', 'rtracklayer'))"

COPY scripts/* /usr/local/bin/
