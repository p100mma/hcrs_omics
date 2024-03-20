FROM rocker/r-ver:4.3.1


# Install Conda
RUN apt-get update && apt-get install -y wget 
RUN mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3


RUN ~/miniconda3/bin/conda init bash
RUN ~/miniconda3/bin/conda init zsh
SHELL ["/bin/bash", "-i", "-c"]
RUN conda install -c bioconda mcl
# Install MCL via Conda
#RUN conda install -c bioconda mcl

RUN install2.r --error --skipinstalled --ncpus -1 \
    Hmisc \
    data.table \
    dynamicTreeCut \
    rmetalog 

RUN mkdir /home/hcrs_omics
#RUN mkdir /home/MCL


#RUN wget -O /home/MCL/soft.gz  'https://micans.org/mcl/src/mcl-14-137.tar.gz' \
#&& tar -xf '/home/MCL/soft.gz' -C /home/MCL 

WORKDIR /home/hcrs_omics
CMD ["/bin/bash"]
