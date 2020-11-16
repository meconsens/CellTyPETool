
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellTyPETool

<!-- badges: start -->

<!-- badges: end -->

\#\#Description The goal of CellTyPETool is to facilitate research that
requires deconvolution of bulk-tissue RNAseq data. The package can be
used to generate cell type proportion estimations from bulk-tissue
RNAseq using two different validated methods, markerGeneProfile and
BRETIGEA. While other R packages have been published that calculate cell
type proportions (McKenzie et al, 2018; Mancarci et al., 2017), they
only offer their own analysis and not a comparison between, nor do they
offer visualizations between cell-type proportions associations with
disease or other important phenotypes. Additionally, they don’t offer
simple way to check which marker genes were encorporated in the
cell-type proportion calculation.

\#\#Installation

To install the latest version of the package:

require(“devtools”) devtools::install\_github(“meconsens/CellTyPETool”,
build\_vignettes = TRUE) library(“CellTyPETool”)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CellTyPETool)
## basic example code
```

\#\#Acknowledgements This package was developed as part of an assessment
for 2020BCB410H: Applied Bioinformatics, University of Toronto,
Toronto,CANADA.

<img src="man/figures/README-pressure-1.png" width="100%" />
