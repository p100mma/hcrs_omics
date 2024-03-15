FROM rocker/r-ver:4.3.1

COPY *.R root/

RUN install2.r --error --skipinstalled --ncpus -1 \
    Hmisc \
    data.table \
    dynamicTreeCut \
    rmetalog 

RUN mkdir /home/hcrs_omics
WORKDIR /home/hcrs_omics
CMD ["/bin/bash"]
