FROM rocker/r-ver:4.3.1

COPY *.R root/

RUN install2.r --error --skipinstalled --ncpus -1 \
    Hmisc \
    data.table \
    dynamicTreeCut \
    rmetalog 

RUN mkdir /home/hcrs_omics
RUN mkdir /home/MCL


RUN wget -O /home/MCL/soft.gz  https://micans.org/mcl/src/mcl-14-137.tar.gz \
&& tar -xf '/home/MCL/soft.gz' -d /home/MCL 

WORKDIR /home/hcrs_omics
CMD ["/bin/bash"]
