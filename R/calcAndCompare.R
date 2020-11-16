#' Calculates and compares the estimations of cell types derived by the
#' markerGeneProfile and BRETIGEA findCells method
#'
#' A function that calculates cell type proportions using two validated methods,
#' markerGeneProfile and findCells, and then compares the cell type proportions
#' calculated by each
#'
#' @param countDf A dataframe with Gene rows and Subject columns.
#' @param bretCellMarkers A dataframe with two columns, markers and cell where
#'  markers are genes for cell types and cell indicates the cell type markers
#'  are the gene for.
#' @param mgpCellMarkers A nested A nested list with as many sub-lists of
#' marker genes as there are for cell types
#'
#' @return Returns a nested list with cell type proportions calculated by the
#' markerGeneProfile method and the BRETIGEA method and a correlation plot to
#' compare the two
#'
#' @examples
#' # Examples 1:
#' # Using countDf, bretCellMarkers, mgpCellMarkers datasets available with package
#'
#' calcAndCompare <- calcAndCompare (
#'                 countDf = countDf,
#'                 mgpCellMarkers = mgpCellMarkers,
#'                 bretCellMarkers = bretCellMarkers)
#'
#' @references
#' H. Wickham. Reshaping data with the reshape package. Journal of Statistical
#' Software, 21(12), 2007.
#'
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
#' York, 2016.
#'
#' Kirill Müller and Hadley Wickham (2020). tibble: Simple Data Frames. R package
#' version 3.0.3. https://CRAN.R-project.org/package=tibble
#'
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific Gene Expression and Coexpression Network Architectures.
#' \emph{Sci Rep} 8, 8868 (2018). https://doi.org/10.1038/s41598-018-27293-5
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
#' URL https://www.R-project.org/.
#'
#' Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
#' https://doi.org/10.21105/joss.01686
#'
#' @export
#' @import markerGeneProfile
#' @import BRETIGEA
#' @import tibble
#' @importFrom dplyr rename
#' @importFrom stats cor
#' @importFrom reshape melt
#' @import ggplot2
#' @importFrom magrittr %>%
calcAndCompare <-function(countDf, mgpCellMarkers, bretCellMarkers){

  if(!all(c("Gene") %in% colnames(countDf))){
    stop("The countDf argument must be a df with a column named Gene (gene symbols)")
  }

  if(!all(c("markers", "cell") %in% colnames(bretCellMarkers))){
    stop("The bretCellMarkers argument must be a df with a column named markers(gene symbols) and a column named cell (corresponding cell types).")
  }

  mgpEstimations<- markerGeneProfile::mgpEstimate(exprData=countDf,
                                                  genes=mgpCellMarkers,
                                                  geneColName="Gene",
                                                  outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                                  geneTransform =NULL, #function(x){homologene::mouse2human(x)$humanGene}, # this is the default option for geneTransform
                                                  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                                                  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                                  removeMinority = TRUE)
  mgp <-mgpEstimations$estimates %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]
  bretEstimations<- BRETIGEA::findCells(countDf, bretCellMarkers, nMarker = 1000, method = "SVD",
                                        scale = TRUE)
  bret <- bretEstimations %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
  mgpFinal <- mgp
  bretFinal <- bret
  colnames(mgp) <- paste("MGP", colnames(mgp), sep = "_")
  mgp<- mgp %>%
    dplyr::rename(
      Sample = MGP_Sample,
    )
  colnames(bret) <- paste("BRET", colnames(bret), sep = "_")
  bret<- bret %>%
    dplyr::rename(
      Sample = BRET_Sample,
    )
  compareEstimates <- merge(mgp, bret, by="Sample")
  compareEstimates <- compareEstimates[,-1]
  markerMat <- data.matrix(compareEstimates, rownames.force = NA)
  cormat <- round(stats::cor(markerMat),2)
  meltedCormat <- reshape::melt(cormat)
  names(meltedCormat)[3] <- "correlation"
  corrPlot <- ggplot2::ggplot(data = meltedCormat, ggplot2::aes(x=X1, y=X2, fill=correlation)) +
    ggplot2::geom_tile() + ggplot2::theme(axis.title.y=element_blank(),
                                          axis.title.x=element_blank(),
                                          axis.text.x = element_text(angle = 45, hjust = 1)) + ggplot2::geom_text(data=meltedCormat,aes(x=X1, y=X2,label=correlation)) +  ggplot2::scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar")

  return(list("mgp" = mgpFinal, "bret" = bretFinal, "corrPlot"= corrPlot))
}
