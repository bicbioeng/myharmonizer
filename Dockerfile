FROM rocker/r-ubuntu:20.04

WORKDIR /app

USER root

ENV PATH=/root/miniconda3/envs/myharmonizer/bin:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update

RUN apt-get install -y wget

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
# update indices
RUN apt-get update -qq
# install two helper packages we need
RUN apt-get install -y --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN pip install ipykernel

COPY myharmonizer .

RUN conda env create -f myharmonizer.yml python=3.10

RUN /root/miniconda3/bin/python -m ipykernel install --user --name myharmonizer

# RUN apt-get install -y libssl-dev

# RUN apt-get install -y libcurl4-openssl-dev

# RUN Rscript -e "install.packages('BiocManager', repos='http://cran.us.r-project.org')"
# RUN Rscript -e "BiocManager::install('RCurl')"
# RUN Rscript -e "BiocManager::install('edgeR')"
# RUN Rscript -e "BiocManager::install('DESeq2')"
# RUN Rscript -e "install.packages('argparse',repos='http://cran.us.r-project.org')"
# RUN Rscript -e "install.packages('plyr',repos='http://cran.us.r-project.org', dependencies = TRUE)"

ENTRYPOINT ["tail", "-f", "/dev/null"]