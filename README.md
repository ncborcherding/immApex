# Apex
API for single-cell immune repertoire deep learning models

<img align="right" src="https://github.com/ncborcherding/Apex/blob/main/www/apex_hex.png" width="305" height="352">

<!-- badges: start -->
[![R-CMD-check](https://github.com/ncborcherding/Apex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/Apex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/Apex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/Apex?branch=main)
<!-- badges: end -->

## Introduction
Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/ncborcherding/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. Further we developed models for embedding TCR sequences ([Trex](https://github.com/ncborcherding/Trex)) 
and BCR sequences ([Ibex](https://github.com/ncborcherding/Ibex)) using convolutional neural networks. Apex is the API for preparing the sequence data for the current and future models
in the scRepertoire ecosystem. 

# System requirements 

Apex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have 
had BCR/TCRs added using [scRepertoire](https://github.com/ncborcherding/scRepertoire). Apex has been tested on OS X and Linux platforms.

**keras** is necessary for Apex (this includes the set up of the tensorflow environment in R):

```r
##Install keras
install.packages("keras")

##Setting up Tensor Flow
library(reticulate)
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(tensorflow)
install_tensorflow()
```

***
### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 
