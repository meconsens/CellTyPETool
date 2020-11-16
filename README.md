
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

``` r
require("devtools")
devtools::install_github("meconsens/CellTyPETool", build_vignettes = TRUE)
library("CellTyPETool")
```

\#\#Overview

``` r
ls("package:CellTyPETool")
```

CellTyPETool contains 4 functions.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CellTyPETool)
#> Warning: replacing previous import 'dplyr::rename' by 'reshape::rename' when
#> loading 'CellTyPETool'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'CellTyPETool'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'CellTyPETool'
#> Warning: replacing previous import 'reshape::expand' by 'tidyr::expand' when
#> loading 'CellTyPETool'
#> Warning: replacing previous import 'magrittr::extract' by 'tidyr::extract' when
#> loading 'CellTyPETool'
## basic example code
```

\#\#Acknowledgements This package was developed as part of an assessment
for 2020 BCB410H: Applied Bioinformatics, University of Toronto,
Toronto,CANADA.

<img src="man/figures/README-pressure-1.png" width="100%" />
