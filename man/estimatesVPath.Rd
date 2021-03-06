% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimatesVPath.R
\name{estimatesVPath}
\alias{estimatesVPath}
\title{Displays strength of association between pathologies and cell types}
\usage{
estimatesVPath(estimates, metadata, cellTypeNames, covar, pathologyName)
}
\arguments{
\item{estimates}{The estimates of cell type proportions returned by calcAndCompare()
either markerGeneProfile derived or BRETIGEA derived}

\item{metadata}{A dataframe with subjects also in countDf and rows
indicating the subjects id, some covariate, and disease state score or
pathology.}

\item{cellTypeNames}{The names of all the unique cell types for which
there are markers in bretCellMarkers: unique(bretCellMarkers$cell)}

\item{covar}{A covariate to be taken into account when running linear models
to check the association between the cell type indicated by cell and the
pathology indicated by pathologyName.}

\item{pathologyName}{The pathology associated with the disease in question
for which the association between it and the cell type indicated by cell is
being examined.}
}
\value{
Returns a volcano plot showing the significance of associations between
each cell type proportion derived and the pathology in question
}
\description{
A function that generates a volcano plot showing the significance of associations between
each cell type proportion derived and the pathology in question
}
\examples{
# Examples 1:
# Using countDf, bretCellMarkers, mgpCellMarkers datasets available with package

calcAndCompareResults <- calcAndCompare (
                countDf = countDf,
                mgpCellMarkers = mgpCellMarkers,
                bretCellMarkers = bretCellMarkers)

estimatesVPathResults <- estimatesVPath(
                estimates = calcAndCompareResults$bret,
                metadata = metadata,
                cellTypeNames = unique(bretCellMarkers$cell),
                covar = "Covariate",
                pathology = "DiseasePhenotypeScore")

}
\references{
H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
York, 2016.

Kamil Slowikowski (2020). ggrepel: Automatically Position Non-Overlapping Text
Labels with 'ggplot2'. R package version 0.8.2.
https://CRAN.R-project.org/package=ggrepel

Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. https://doi.org/10.1523/ENEURO.0212-17.201

McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific Gene Expression and Coexpression Network Architectures.
\emph{Sci Rep} 8, 8868 (2018). https://doi.org/10.1038/s41598-018-27293-5

R Core Team (2020). R: A language and environment for statistical computing. R
Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
https://doi.org/10.21105/joss.01686

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
}
