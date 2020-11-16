
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellTyPETool

<!-- badges: start -->

<!-- badges: end -->

\#\#Description

The goal of CellTyPETool is to facilitate research that requires
deconvolution of bulk-tissue RNAseq data. The package can be used to
generate cell type proportion estimations from bulk-tissue RNAseq using
two different validated methods, markerGeneProfile and BRETIGEA. While
other R packages have been published that calculate cell type
proportions (McKenzie et al, 2018; Mancarci et al., 2017), they only
offer their own analysis and not a comparison between, nor do they offer
visualizations between cell-type proportions associations with disease
or other important phenotypes. Additionally, they don’t offer simple way
to check which marker genes were encorporated in the cell-type
proportion calculation.

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
## basic example code
```

\#\#References

H. Wickham. Reshaping data with the reshape package. Journal of
Statistical Software, 21(12), 2007.

Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for
Visualization. R package version 1.1.1.
<https://CRAN.R-project.org/package=scales>

Kamil Slowikowski (2020). ggrepel: Automatically Position
Non-Overlapping Text Labels with ‘ggplot2’.R package version 0.8.2.
<https://CRAN.R-project.org/package=ggrepel>

Kirill Müller and Hadley Wickham (2020). tibble: Simple Data Frames. R
package version 3.0.3. <https://CRAN.R-project.org/package=tibble>

Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille,
E., & Pavlidis, P. (2017). CrossLaboratory Analysis of Brain Cell Type
Transcriptomes with Applications to Interpretation of Bulk Tissue Data.
eNeuro, 4(6), ENEURO.0212-17.2017.
<https://doi.org/10.1523/ENEURO.0212-17.2017>

McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression
analysis of multifactor RNA-Seq experiments with respect to biological
variation. Nucleic Acids Research 40, 4288-4297

McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific
Gene Expression and Coexpression Network Architectures. Sci Rep 8, 8868
(2018). <https://doi.org/10.1038/s41598-018-27293-5>

R Core Team (2020). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. URL
<https://www.R-project.org/>.

Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe
Operator for R. R package version 1.5.
<https://CRAN.R-project.org/package=magrittr>

Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source
Software, 4(43), 1686, <https://doi.org/10.21105/joss.01686>

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York. ISBN 978-3-319-24277-4,
<https://ggplot2.tidyverse.org>.

\#\#Acknowledgements This package was developed as part of an assessment
for 2020 BCB410H: Applied Bioinformatics, University of Toronto,
Toronto,CANADA.

<img src="man/figures/README-pressure-1.png" width="100%" />
