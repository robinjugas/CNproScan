FROM bioconductor/bioconductor_docker:RELEASE_3_13

LABEL maintainer "Robin Jugas <robinjugas@gmail.com>"
LABEL Description="CNproScan The CNV detection R package" Vendor="BUT BRNO 2020" Version="1.0"

# INSTALL PACKAGES

RUN R -e "install.packages(c('parallel','foreach','doParallel','xlsx', 'devtools'))"

RUN R -e "BiocManager::install(c('seqinr','Rsamtools','GenomicRanges','IRanges'))"

RUN R -e "devtools::install_github('robinjugas/CNproScan')"
