#' Returns genes used, removedMarkerRatios and proportion of variance explained
#' by the first PC per cell type in the Marker Gene Profile method
#'
#' A function that returns a df of the genes used, removedMarkerRatios and proportion
#' of variance explained by the first PC in the markerGeneProfile estimation method
#'
#' @param countDf A dataframe with a Gene row in HGNC symbols and then Subject columns.
#' @param mgpCellMarkers A nested list with as many sub-lists of
#' marker genes as there are for cell types
#'
#' @return Returns a dataframe with the cell type, markersUsed (list of marker genes used
#' per cell type), removedMarkerRatios (list of removed marker
#' ratios per cell type) and percentVariancePC1 (list of variance explained
#' by the first PC per cell type)
#'
#'@examples
#' # Examples 1:
#' # Using countDf and mgpCellMarkers datasets available with package
#'
#' mgpQCMetricsResults <- mgpQCMetrics(countDf, mgpCellMarkers)
#'
#' @references
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017.https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe Operator
#' for R. R package version 1.5. https://CRAN.R-project.org/package=magrittr
#' @export
#'
mgpQCMetrics <-function(countDf, mgpCellMarkers){

  if(!all(c("Gene") %in% colnames(countDf))){
    stop("The countDf argument must be a df with a column named Gene (HGNC gene symbols)")
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
    masterlist <- paste0(rownames(cellsDf), collapse=', ')
    rmMarkerRatios <- mgpEstimations$removedMarkerRatios[i]
    percentVariance <- ((summary(mgpEstimations$trimmedPCAs[[i]]))[6]) %>% as.data.frame()
    percentVariancePC1 <- percentVariance[2,1]
    if(i==1){
      masterDf <- data.frame( "markersUsed" = masterlist, "removedMarkerRatios" = rmMarkerRatios,
                              "percentVariancePC1" = percentVariancePC1)
    }
    else{
      df <- data.frame( "markersUsed" = masterlist, "removedMarkerRatios" = rmMarkerRatios,
                        "percentVariancePC1" = percentVariancePC1)
      masterDf <- rbind(masterDf, df)
    }
  }
  masterDf <- tibble::rownames_to_column(masterDf, var = "cellType")
  return(masterDf)
}
