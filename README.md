# myharmonizer: Python package to harmonize and assess similarity between new user data and an existing knowledge base in a myHarmonizer object

In order to facilitate the comparison between an existing bulk RNA-seq knowledge base and new user data, the myharmonizer package was designed to take
a myHarmonizer object output from the [DeepSeqDock framework](https://github.com/bicbioeng/DeepSeqDock) and use the frozen preprocessing, 
scaling, and transformation methods used to build the knowledge base on new user data. By applying these methods on new user data, the new data 
is brought into a similar data representation as the existing knowledge base and similarity can be calculated within a theoretically more meaningful
space. 

In short, myHarmonizer performs the following:
* Transform new user datasets into the same representation as the input myHarmonizer knowledge base
* Calculate similarities between user samples and knowledge samples
* Visualize similarity matrices

## Installation
### Using Docker

The most straightforward way of running myharmonizer via command line is by setting up a Docker instance. First, make sure Docker is installed, then run:

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
 
### Install package directly in Ubuntu

Requirements:
 - python <=3.9.12  

Clone git repo
```
git clone https://github.com/bicbioeng/myharmonizer
```
-or- install package using pip
```
pip install ./myharmonizer
```

Install R in Unbuntu (if not previously installed)
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

## Running myharmonizer package
### Import packages and toy myHarmonizer object
```
python

import pandas as pd
import numpy as np
import myharmonizer as mh
from importlib.resources import files


# Build toy myHarmonizer
toy_mH = mh.myHarmonizer(files(myharmonizer) / 'myHarmonizer-YYYY_mm_dd-HH_MM_SS.json')
```

### Generate synthetic data
```
# Get the raw feature list for the toy dataset
rawfeatures = toy_mH.modelmeta['encoder_metrics']['features']

# Generate 5 random samples
newdata = pd.DataFrame(np.random.poisson(size=(5,len(rawfeatures))), 
    index=['A', 'B', 'C', 'D', 'E'],
    columns=rawfeatures)
```

### Transform synthetic data to representation space of toy myHarmonizer knowledge base
```
transformeddata = toy_mH.transform(newdata)
```
This step can also be broken down into individual transformations for niceify (harmonize feature list,
outlier removal, and add one psuedocount), normalize, scale, and encode.

### Calculate similarities

The representation of the knowledge base is stored as the data attribute of the myHarmonizer class. Consequently,
to compare the new data to the existing data representation, the following can be used:

```
pearson_sim = mh.similarity(transformeddata, toy_mH.data, metric='Pearson')
```

For this function, the similarity metric can be selected from one of 'Pearson', 'Spearman', 'CCC', 'Euclidean', 'Manhattan', 
or 'Cosine'.

### Visualize similarities

To visualize the similarities as a heatmap (scaled from 0-1), there is a convenience function that implements a heatmaps based on the
seaborn clustermap function. Users can imput their own metadata for their sample index as a Series, and/or select one of the columns 
of the metadata in the knowledge base, if metadata is included in the myHarmonizer object.

```
# Examine metadata in myHarmonizer object
toy_mH.metadata

# Plot heatmap
mh.heatmap(pearson_sim, toy_mH, user_metadata=None, kb_metadata='Meta A')
```


