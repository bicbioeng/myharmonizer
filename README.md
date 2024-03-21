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

These functions can also be performed in the online application at: [https://myharmonizer.bicbioeng.org](https://myharmonizer.bicbioeng.org)

## Installation
### Using Docker

The most straightforward way of running myharmonizer via command line is by setting up a Docker instance. First, make sure Docker is installed. If using a Windows system, then the Docker can run through WSL2. Next, download the Docker image from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10119419.svg)](https://doi.org/10.5281/zenodo.10119419)

Load the docker image:

```
docker load < myharmonizerdock.tar.gz
```

then run the container interactively:

```
docker run -it --rm -v "$HOME"/myharmonizeroutput:/app/myharmonizeroutput myharmonizerdock
```
Parameters:
 - -it: interactive mode
 - --rm: automatically clean container when container stopped
 - -v: bind container volume to host folder. Host folder will by default be in myharmonizeroutput folder of WSL2 or Linux home folder.
 
 Scripts are not automatically added to the path, so it is easiest to run python in the default /app/myharmonizer directory. 
 
### Install package with conda environment

Requirements:
 - python <=3.9.12  

Clone git repo
```
git clone https://github.com/bicbioeng/myharmonizer
```
In the myharmonizer/myharmonizer directory, create a conda environment from the .yml file.

```
conda env create -f myharmonizer.yml
```

and then activate the environment before running the modules.

```
conda activate myharmonizer
```

For Ubuntu version 20.10 and above, it may be necessary to install libffi7_3. The dependency chain with R does not allow for updated versions of python to be used at the time of this writing.

## Running myharmonizer package
### 1) Import packages and toy myHarmonizer object

In python:
```python
import pandas as pd
import numpy as np
import myharmonizer as mh
from matplotlib import pyplot as plt


# Build toy myHarmonizer
toy_mH = mh.myHarmonizer('myHarmonizer-YYYY_mm_dd-HH_MM_SS.json')
```

### 2) Generate synthetic data
```python
# Get a shortened feature list for the toy dataset
rawfeatures = toy_mH.modelmeta['encoder_metrics']['features']

# Generate 5 random samples
newdata = pd.DataFrame(np.random.poisson(size=(5,len(rawfeatures))), 
    index=['A', 'B', 'C', 'D', 'E'],
    columns=rawfeatures)
```

When running read user data, please format data as a pandas DataFrame with samples as rows
and unique features as columns. While the feature list will be reconciled, users should choose a 
feature identifier (e.g. ENSEMBL, Official Gene Symbol) that most closely aligns with that
of the knowledge base. The original feature list for the knowledge base can be conveniently
accessed as the 'features' attribute of the myHarmonizer class if unknown:

```python
toy_mH.features
```

### 3) Transform synthetic data to representation space of toy myHarmonizer knowledge base
```python
transformeddata = toy_mH.transform(newdata)
```
This step can also be broken down into individual transformations for niceify (harmonize feature list,
outlier removal, and add one psuedocount), normalize, scale, and encode.

### 4) Calculate similarities

The representation of the knowledge base is stored as the data attribute of the myHarmonizer class. Consequently,
to compare the new data to the existing data representation, the following can be used:

```python
pearson_sim = mh.similarity(transformeddata, toy_mH.data, metric='Pearson')
```

For this function, the similarity metric can be selected from one of 'Pearson', 'Spearman', 'CCC', 'Euclidean', 'Manhattan', 
or 'Cosine'.

### 5) Visualize similarities

To visualize the similarities as a heatmap (scaled from 0-1), there is a convenience function that implements a heatmaps based on the
seaborn clustermap function. Users can imput their own metadata for their sample index as a Series, and/or select one of the columns 
of the metadata in the knowledge base, if metadata is included in the myHarmonizer object.

```python
# Examine metadata in myHarmonizer object
toy_mH.metadata

# Plot heatmap
mh.heatmap(pearson_sim, toy_mH, user_metadata=None, kb_metadata='Meta A')
plt.savefig('Fig.png')
```

## Citations and Licensing

myharmonizer: a Python package for the transformation of add-on datasets and similarity calculations with datasets in myharmonizer JSON objects. \
Copyright (C) 2024 bicbioeng

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

This work includes functions modified from the calcNormFactors function in edgeR and the estimateSizeFactors function from DESeq2. Modified functions are marked and can be found in the GeTMM_preprocessing.R and GeVST_preprocessing.R files. Original code can be found under calcNormFactors_edgeR.R ([https://code.bioconductor.org/browse/edgeR/blob/RELEASE_3_16/R/calcNormFactors.R](https://code.bioconductor.org/browse/edgeR/blob/RELEASE_3_16/R/calcNormFactors.R)) or estimateSizeFactorsForMatrix_DESeq2.R ([https://code.bioconductor.org/browse/DESeq2/blob/RELEASE_3_12/R/core.R]( https://code.bioconductor.org/browse/DESeq2/blob/RELEASE_3_12/R/core.R)). DESeq2 is distributed under the [LGPL license (>=3)](https://www.gnu.org/licenses/lgpl-3.0.en.html) and edgeR is distributed under the [LGPL license (>=2)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html). A copy of LGPL 3.0 is also available in this repository. 

DESeq2 Citation:

Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)  

edgeR Citations:

  1) Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
  expression analysis of digital gene expression data. Bioinformatics 26, 139-140

  2) McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq
  experiments with respect to biological variation. Nucleic Acids Research 40, 4288-4297

  3) Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways: differential expression
  analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
  F1000Research 5, 1438

  4) Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2024). edgeR 4.0: powerful differential analysis
  of sequencing data with expanded functionality and improved support for small counts and larger
  datasets. bioRxiv doi: 10.1101/2024.01.21.576131
