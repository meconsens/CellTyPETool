#' Determines associations between cell type estimates and pathology
#'
#' A function that runs a linear model on the cell-type
#' proportion estimates by BRETIGEA.
#'
#' @param markerMat The $SPVS of the return of findCellsMod once converted to a
#'    dataframe with the first column being Sample corresponding to the metadata
#'    Sample id. (Dataframe with Sample column and then estimate of cell type
#'    proportion calculated by BRETIGEA findCells for each cell type in
#'    cellTypeNames)
#' @param metadata A dataframe with subjects also in countDf and rows
#'    indicating the subjects id, some covariate, and disease state score or
#'    pathology.
#' @param covar A covariate to be taken into account when running linear models
#'     to check the association between the cell type indicated by cell and the
#'     pathology indicated by pathologyNames.
#' @param pathologyName The pathology associated with the disease in question
#'     for which the association between it and the cell type indicated by cell is
#'     being examined.
#' @param cellTypeNames The names of all the unique cell types for which
#'     there are markers in bretCellMarkers: unique(bretCellMarkers$cell)
#'
#' @return Returns a matrix array ready to be formatted by bretigeaMarkersAddition
#'
#' @examples
#' # Examples 1:
#' # Using bretCellMarkers, metadata datasets available with package
#'
#' rownames(countDf) <- countDf$Gene
#' countDf <- countDf[,-1]
#' bretigeaMarkers <- findCellsMod(
#'                inputMat = countDf,
#'                markers = bretCellMarkers,
#'                nMarker = 10,
#'                method = "SVD",
#'                scale = TRUE)
#' markerMat <- as.data.frame(bretigeaMarkers$SPVS)
#' markerMat <- tibble::rownames_to_column(markerMat, var = 'Sample')
#' markersPathologyResults <- markersPathology(
#'                       markerMat = markerMat,
#'                       metadata = metadata,
#'                       covar = "Covariate",
#'                       pathologyName = "DiseasePhenotypeScore",
#'                       cellTypeNames = unique(bretCellMarkers$cell))
#'
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#'
#' Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}, Wadsworth & Brooks/Cole.
#'
#' R Core Team (2020). R: A language and environment for statistical computing.
#' R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#'
#' Wilkinson, G. N. and Rogers, C. E. (1973). Symbolic descriptions of factorial models for analysis of variance. \emph{Applied Statistics}, 22, 392–399. doi: 10.2307/2346786.
#'
#' @export
#' @importFrom stats anova coef lm as.formula
#' @importFrom magrittr %>%
markersPathology <-function(markerMat, metadata, covar, pathologyName, cellTypeNames){

  if(!covar %in% colnames(metadata)) stop("The covar argument must have a corresponding column in metadata.")

  if(!pathologyName %in% colnames(metadata)) stop("The pathologyName argument must have a corresponding column in metadata.")

  markersMeta <- merge(markerMat, metadata, by="Sample")
  model.data <- markersMeta
  #get the direction of effect and significance of association between cell type and pathology, as well as number of observations
  results <- sapply(cellTypeNames,function(celltype) {
    sapply(pathologyName, function(pathology) {
      form <- as.formula(paste0(celltype,"~",pathology," + ",paste0(covar,collapse=" + ")))
      model <- stats::lm(data=model.data,form)
      p <- stats::anova(model)[pathology,5]
      beta <- stats::coef(model)[pathology]
      n <- nrow(model$model)
      return(c(celltype,pathology,beta, p,n))
    })
  })
}



#' Runs BRETIGEA using top 1, top 1 & 2, top 1,2,3 ... to top 1...n markers
#' A function that reruns the BRETIGEA findCells method (modified) on a loop to see
#' the influence of different combinations of markers in determining the cell-type
#' proportion estimate.
#'
#' @param countDf A dataframe with Gene rows and Subject columns.
#' @param bretCellMarkers A dataframe with two columns, markers and cell where
#'      markers are genes for cell types and cell indicates the cell type markers
#'      are the gene for.
#' @param cellName A string indicating which of the cells in the cell type marker list
#'     is being specified to look at.
#' @param covar A covariate to be taken into account when running linear models
#'     to check the association between the cell type indicated by cell and the
#'     pathology indicated by pathologyNames.
#' @param metadata A dataframe with subjects also in countDf and rows
#'     indicating the subjects id, some covariate, and disease state score or
#'     pathology.
#' @param pathologyName The pathology associated with the disease in question
#'     for which the association between it and the cell type indicated by cell is
#'     being examined.
#' @param n Specifies how many times BRETIGEA will run BRETIGEA using top 1,
#'     top 1 & 2, top 1,2,3 ... to top 1...n markers
#' @param cellTypeNames The names of all the unique cell types for which
#'     there are markers in bretCellMarkers: unique(bretCellMarkers$cell)
#'
#' @return Returns a graph of the significance of the cell type proportion
#' specified's association to the pathology indicated upon marker addition
#' from 0 to n
#'
#' @examples
#' # Examples 1:
#' # Using countDf, bretCellMarkers, metadata datasets available with package
#'
#' bretMarkerEffectOnPathologyResults <- bretMarkerEffectOnPathology (
#'                           countDf = countDf,
#'                           bretCellMarkers = bretCellMarkers,
#'                           cellName = "mic",
#'                           metadata = metadata,
#'                           covar = "Covariate",
#'                           pathologyName = "DiseasePhenotypeScore",
#'                           cellTypeNames = unique(bretCellMarkers$cell),
#'                           n= 10)
#'
#' @references
#' Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization.
#' R package version 1.1.1. https://CRAN.R-project.org/package=scales
#'
#' Kirill Müller and Hadley Wickham (2020). tibble: Simple Data Frames. R package
#' version 3.0.3. https://CRAN.R-project.org/package=tibble
#'
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#'
#' Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe Operator for R.
#' \emph{R package version 1.5}.https://CRAN.R-project.org/package=magrittr
#'
#' Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
#' https://doi.org/10.21105/joss.01686
#'
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr filter
#' @import tibble
#' @import scales
#' @importFrom stats p.adjust
#' @import utils
#' @import ggplot2
#' @import ggrepel
#' @importFrom magrittr %>%

bretMarkerEffectOnPathology <- function(countDf, bretCellMarkers, cellName, metadata, covar, pathologyName, cellTypeNames, n){
  if(!all(c("markers", "cell") %in% colnames(bretCellMarkers))){
    stop("The bretCellMarkers argument must be a df with a column named marker s(gene symbols) and a column named cell (corresponding cell types).")
  }
  if(!cellName %in% cellTypeNames) stop("The cellName argument must specify a cell in cellTypeNames.")

  if(!covar %in% colnames(metadata)) stop("The covar argument must have a corresponding column in metadata.")

  if(!pathologyName %in% colnames(metadata)) stop("The pathologyName argument must have a corresponding column in metadata.")

  if(n > 1000) {
    warning("n is too large", call. = FALSE)
    n<- 1000 # correct the input for user
  }

  #prep countDf for use by BRETIGEA
  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]
  for (i in 1: n)
  {
    #run the modified BRETIGEA method for calculating cell type proportion estimates, but get the markers used back
    bretigeaMarkers <- findCellsMod(countDf, bretCellMarkers, nMarker = i, method = "SVD",scale = TRUE)
    bretigeaMarkers$SPVS <- bretigeaMarkers$SPVS %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
    #get direction of effect significance of association between cell type and pathology, as well as number of observations
    results <- markersPathology(bretigeaMarkers$SPVS, metadata, covar, pathologyName, cellTypeNames)
    results <- as.data.frame(matrix(results,ncol=5,byrow=T),stringsAsFactors = F)
    names(results) <- c("celltype","pathology","beta","p","n")
    #correctly assign types of results
    results <- within(results,{
      p <- as.numeric(p)
      beta <- as.numeric(beta)
      n <- as.numeric(n)})
    results$fdr <- stats::p.adjust(results$p, method="fdr")
    results$nMarker <- i
    if(i > 1){
      results_final <- rbind(results_final, results)
    }
    else{
      results_final <- results
    }
  }
  #get the information from linear model about specific cell type
  cellPathology <- results_final %>% dplyr::filter(celltype == cellName)
  cellPathology <- cellPathology %>% tidyr::gather(key, value, pathology)
  #get cell type marker gene names
  cellLabel<- bretigeaMarkers$markerList%>% dplyr::filter(cell == cellName)
  cellLabel <- utils::head(cellLabel,n)
  cellPathology$markers <- cellLabel$markers
  #make the line plot, labelling marker added at each observation
  graph <- cellPathology %>%
    ggplot2::ggplot(ggplot2::aes(x=nMarker, y=-log10(p), color=value)) + ggplot2::theme_minimal() + ggrepel::geom_text_repel(aes(label=markers), size =2)+ ggplot2::geom_line(aes(color = value), size =1) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      colour = 'black', size = 13,
      hjust = 0.5, vjust = 0.5)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
      colour = 'black', size = 13,
      hjust = 0.5, vjust = 0.5)) + ggplot2::ggtitle(paste0("Significance of ", cellName,  " Cell Type Proportion Association to ", pathologyName, " Upon Marker Addition"))
  return(graph)
}

#[END]
