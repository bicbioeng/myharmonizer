# myharmonizer
Python package to harmonize and evaluate similarity between new data and an existing knowledge base in a myHarmonizer object.

## Installation
### Using Docker
```
docker run -d -e GRANT_SUDO=yes --user root --rm -p 8888:8888 -e NB_USER='jovyan' -e CHOWN_HOME=yes -w "/app" -e CHOWN_HOME_OPTS='-R' -v /home/tuyendo/USD/app/:/app --name Myharmonizer us-central1-docker.pkg.dev/nosi-usd-biofilm/nosi-usd-biofilm-arti/myharmonizer
```
Parameters:
 - -d: detach mode
 - GRANT_SUDO: grant super user for nb user
 - --user: specify docker execute user
 - --rm: automatic clean container when container stop
 - -p: export container container port for host using [host:container]
 - NB_USER: specify nb user
 - CHOWN_HOME=yes -w "/app" -e CHOWN_HOME_OPTS='-R' : Change working directory
 - -v : bind container volume to host folder
 - --name : specify container name

Requirements:
 - python <=3.9.12  

Clone git repo
```
git clone https://github.com/bicbioeng/myharmonizer
```
Install package using pip
```
pip install ./myharmonizer
```
### Install R dependencies in Ubuntu
Install R in Unbuntu
```
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
```
Install R dependencies in Ubuntu
```
apt-get install -y libssl-dev
apt-get install -y libcurl4-openssl-dev

Rscript -e "install.packages('BiocManager', repos='http://cran.us.r-project.org')"
Rscript -e "BiocManager::install('RCurl')"
Rscript -e "BiocManager::install('edgeR')"
Rscript -e "BiocManager::install('DESeq2')"
Rscript -e "install.packages('argparse',repos='http://cran.us.r-project.org')"
Rscript -e "install.packages('plyr',repos='http://cran.us.r-project.org', dependencies = TRUE)"
```