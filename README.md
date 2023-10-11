# scGOclust
## Leveraging Gene Ontology to Measure Cell Type Similarity Between Single Cell RNA-Seq Datasets

Author & Maintainer: Yuyao Song <ysong@ebi.ac.uk>

![cran downloads](https://cranlogs.r-pkg.org/badges/grand-total/scGOclust)

### Installation

1. First create the conda environment. Mamba is recommended as a faster alternative for conda.

    `conda env create -f scGOclust_conda_7Dec2022.yml`

2. Then, open R under this environment, and install several packages not in conda: 

    `remotes::install_github('satijalab/seurat-wrappers'), install.packages("pheatmap", "slanter")`


3. Finally, install `scGOclust` from GitHub: 
    
    `devtools::install_github("YY-SONG0718/scGOclust", ref = "main")`

### Usage

This package operates on pairs of `Seurat` objects

Refer to `vignettes` for usage examples


### License: GPLv3.0
