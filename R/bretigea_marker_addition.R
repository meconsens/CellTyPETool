#' Determines associations between cell type estimates and pathology
#'
#' A function that runs a linear model on the cell-type
#' proportion estimates by BRETIGEA.
#'
#' @param markerMat The $SPVS of the return of find_cells_mod once converted to a
#' dataframe with the first column being Sample corresponding to the metadata
#' Sample id. (Dataframe with Sample column and then estimate of cell type
#' proportion calculated by BRETIGEA findCells for each cell type in
#' cell_type_names)
#' @param metadata A dataframe with subjects also in count_df and rows
#' indicating the subjects id, some covariate, and disease state score or
#' pathology.
#' @param covar A covariate to be taken into account when running linear models
#' to check the association between the cell type indicated by cell and the
#' pathology indicated by pathology_names.
#' @param pathology_name The pathology associated with the disease in question
#' for which the association between it and the cell type indicated by cell is
#' being examined.
#' @param cell_type_names The names of all the unique cell types for which
#' there are markers in bret_cell_markers: unique(bret_cell_markers$cell)
#'
#' @return Returns a matrix array ready to be formatted by
#'
#' @examples
#' # Examples 1:
#' # Using bret_cell_markers, metadata datasets available with package
#'
#' rownames(count_df) <- count_df$Gene
#' count_df <- count_df[,-1]
#' bretigeaMarkers <- find_cells_mod(
#'                inputMat = count_df,
#'                markers = bret_cell_markers,
#'                nMarker = 10,
#'                method = "SVD",
#'                scale = TRUE)
#' markerMat <- as.data.frame(bretigeaMarkers$SPVS)
#' markerMat <- tibble::rownames_to_column(markerMat, var = 'Sample')
#' markers_pathology <- markers_pathology(
#'                       markerMat = markerMat,
#'                       metadata = metadata,
#'                       covar = "Covariate",
#'                       pathology_name = "DiseasePhenotypeScore",
#'                       cell_type_names = unique(bret_cell_markers$cell))
#'
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#'
#'Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}, Wadsworth & Brooks/Cole.
#'
#'Wilkinson, G. N. and Rogers, C. E. (1973). Symbolic descriptions of factorial models for analysis of variance. \emph{Applied Statistics}, 22, 392–399. doi: 10.2307/2346786.
#'
#'@export
#'@import stats
#'@importFrom magrittr %>%
markers_pathology <-function(markerMat, metadata, covar, pathology_name, cell_type_names){
  markers_meta <- merge(markerMat, metadata, by="Sample")
  model.data <- markers_meta
  results <- sapply(cell_type_names,function(celltype) {
    sapply(pathology_name, function(pathology) {
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
#'
#' A function that reruns the BRETIGEA findCells method (modified) on a loop to see
#' the influence of different combinations of markers in determining the cell-type
#' proportion estimate.
#'
#' @param count_df A dataframe with Gene rows and Subject columns.
#' @param bret_cell_markers A dataframe with two columns, markers and cell where
#'  markers are genes for cell types and cell indicates the cell type markers
#'  are the gene for.
#' @param cell A string indicating which of the cells in the cell type marker list
#' is being specified to look at.
#' @param covar A covariate to be taken into account when running linear models
#' to check the association between the cell type indicated by cell and the
#' pathology indicated by pathology_names.
#' @param metadata A dataframe with subjects also in count_df and rows
#' indicating the subjects id, some covariate, and disease state score or
#' pathology.
#' @param pathology_name The pathology associated with the disease in question
#' for which the association between it and the cell type indicated by cell is
#' being examined.
#' @param n Specifies how many times BRETIGEA will run BRETIGEA using top 1,
#' top 1 & 2, top 1,2,3 ... to top 1...n markers
#' @param cell_type_names The names of all the unique cell types for which
#' there are markers in bret_cell_markers: unique(bret_cell_markers$cell)
#'
#' @return Returns a graph of the significance of the cell type proportion
#' specified's association to the pathology indicated upon marker addition
#' from 0 to n
#'
#' @examples
#' # Examples 1:
#' # Using count_df, bret_cell_markers, metadata datasets available with package
#'
#' bretigea_marker_addition <- bretigea_marker_addition (
#'                           count_df = count_df,
#'                           bret_cell_markers = bret_cell_markers,
#'                           cell = "mic",
#'                           metadata = metadata,
#'                           covar = "Covariate",
#'                           pathology_name = "DiseasePhenotypeScore",
#'                           cell_type_names = unique(bret_cell_markers$cell),
#'                           n= 10)
#'
#' @references
#'Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#'CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#'Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. \href{https://doi.org/10.1523/ENEURO.0212-17.201}
#'
#'Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe Operator for R.
#'\emph{R package version 1.5}. \href{https://CRAN.R-project.org/package=magrittr}
#'
#'Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
#'https://doi.org/10.21105/joss.01686
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import scales
#' @import stats
#' @import utils
#' @import ggplot2
#' @import ggrepel
#' @importFrom magrittr %>%

bretigea_marker_addition <- function(count_df, bret_cell_markers, cell, metadata, covar, pathology_name, cell_type_names, n){
  rownames(count_df) <- count_df$Gene
  count_df <- count_df[,-1]
  for (i in 1: n)
  {
    bretigeaMarkers <- find_cells_mod(count_df, bret_cell_markers, nMarker = i, method = "SVD",scale = TRUE)
    bretigeaMarkers$SPVS <- bretigeaMarkers$SPVS %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
    results <- markers_pathology(bretigeaMarkers$SPVS, metadata, covar, pathology_name, cell_type_names)
    results <- as.data.frame(matrix(results,ncol=5,byrow=T),stringsAsFactors = F)
    names(results) <- c("celltype","pathology","beta","p","n")
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
  cell_pathology <- results_final %>% dplyr::filter(celltype == cell)
  cell_pathology <- cell_pathology %>% tidyr::gather(key, value, pathology)
  cell_label<- bretigeaMarkers$marker_list%>% dplyr::filter(cell == cell)
  cell_label <- utils::head(cell_label,n)
  cell_pathology$markers <- cell_label$markers
  graph <- cell_pathology %>%
    ggplot2::ggplot(ggplot2::aes(x=nMarker, y=-log10(p), color=value)) + ggplot2::theme_minimal() + ggrepel::geom_text_repel(aes(label=markers), size =2)+ ggplot2::geom_line(aes(color = value), size =1) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      colour = 'black', size = 13,
      hjust = 0.5, vjust = 0.5)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
      colour = 'black', size = 13,
      hjust = 0.5, vjust = 0.5)) + ggplot2::ggtitle(paste0("Significance of ", cell,  " Cell Type Proportion Association to ", pathology_name, " Upon Marker Addition"))
  return(graph)
}

