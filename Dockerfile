FROM rocker/r-ver:4.3.1

COPY *.R root/

RUN install2.r --error --skipinstalled --ncpus -1 \
    Hmisc \
    data.table \
    dynamicTreeCut \
    rmetalog 

CMD ["/bin/bash"]
CMD ["cd /home/"]
CMD ["mkdir hcrs_omics"]
CMD ["cd hcrs_omics"]
