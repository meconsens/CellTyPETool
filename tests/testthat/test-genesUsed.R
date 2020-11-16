#' Returns genes used in Marker Gene Profile and BRETIGEA method
#'
#' A function that returns the genes used in the markerGeneProfile estimation method
#' and the BRETIGEA method
#'
#' @param countDf A dataframe with Gene rows and Subject columns.
#' @param bretCellMarkers A dataframe with two columns, markers and cell where
#'  markers are genes for cell types and cell indicates the cell type markers
#'  are the gene for.
#' @param mgpCellMarkers A nested A nested list with as many sub-lists of
#' marker genes as there are for cell types
#'
#' @return Returns a nested list with mgpMarkers and bretMarkers which are each
#' lists of markers used for each method
#'
#' @examples
#' # Examples 1:
#' # Using countDf, bretCellMarkers, mgpCellMarkers datasets available with package
#'
#' genesUsed <- genesUsed (
#'                 countDf = countDf,
#'                 mgpCellMarkers = mgpCellMarkers,
#'                 bretCellMarkers = bretCellMarkers)
#'
#' @references
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017.https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific Gene Expression and Coexpression Network Architectures.
#' \emph{Sci Rep} 8, 8868 (2018). https://doi.org/10.1038/s41598-018-27293-5
#'
#' Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe Operator
#' for R. R package version 1.5. https://CRAN.R-project.org/package=magrittr
#'
#' @export
#' @import markerGeneProfile
#' @importFrom magrittr %>%
genesUsed<-function(countDf, mgpCellMarkers, bretCellMarkers){
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
  i= 0
  for(cell in names(mgpCellMarkers)){
    i = i + 1
    cellsDf <- mgpEstimations$usedMarkerExpression[i] %>% as.data.frame()
    if(i==1){
      masterlist <- data.frame("MarkersUsed" = rownames(cellsDf), "Cell" = cell)
    }
    else{
      list <- data.frame("MarkersUsed" = rownames(cellsDf), "Cell" = cell)
      masterlist <- rbind(masterlist, list)
    }
  }
  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]
  bretEstimations<- findCellsMod(countDf, bretCellMarkers, nMarker = 1000, method = "SVD",
                                 scale = TRUE)
  colnames(bretEstimations$markerList) <- c("MarkersUsed", "Cell")
  return(list("mgpMarkers"= masterlist, "bretMarkers"= bretEstimations$markerList))
}
