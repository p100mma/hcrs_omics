FROM rocker/r-ver:4.3.1


# Install Conda
RUN apt-get update && apt-get install -y wget bzip2 ca-certificates curl git && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Install MCL via Conda
RUN conda install -c bioconda mcl

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
