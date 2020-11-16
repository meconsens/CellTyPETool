#' Calculates and compares the estimations of cell types derived by the
#' markerGeneProfile and BRETIGEA findCells method
#'
#' A function that calculates cell type proportions using two validated methods,
#' markerGeneProfile and findCells, and then compares the cell type proportions
#' calculated by each
#'
#' @param count_df A dataframe with Gene rows and Subject columns.
#' @param bret_cell_markers A dataframe with two columns, markers and cell where
#'  markers are genes for cell types and cell indicates the cell type markers
#'  are the gene for.
#' @param mgp_cell_markers A nested A nested list with as many sub-lists of
#' marker genes as there are for cell types
#'
#' @return Returns a nested list with cell type proportions calculated by the
#' markerGeneProfile method and the BRETIGEA method and a correlation plot to
#' compare the two
#'
#' @examples
#' # Examples 1:
#' # Using count_df, bret_cell_markers, mgp_cell_markers datasets available with package
#'
#' calc_and_compare <- calc_and_compare (
#'                 count_df = count_df,
#'                 mgp_cell_markers = mgp_cell_markers,
#'                 bret_cell_markers = bret_cell_markers)
#'
#' @references
#'Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#'CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#'Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. \href{https://doi.org/10.1523/ENEURO.0212-17.201}
#'
#'McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific Gene Expression and Coexpression Network Architectures.
#'\emph{Sci Rep} 8, 8868 (2018). https://doi.org/10.1038/s41598-018-27293-5
#'
#'#'Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
#'https://doi.org/10.21105/joss.01686
#'
#' @export
#' @import markerGeneProfile
#' @import BRETIGEA
#' @import tibble
#' @import dplyr
#' @import stats
#' @import reshape
#' @import ggplot2
calc_and_compare <-function(count_df, mgp_cell_markers, bret_cell_markers){
  mgp_estimations<- markerGeneProfile::mgpEstimate(exprData=count_df,
                                                   genes=mgp_cell_markers,
                                                   geneColName="Gene",
                                                   outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                                   geneTransform =NULL, #function(x){homologene::mouse2human(x)$humanGene}, # this is the default option for geneTransform
                                                   groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                                                   seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                                   removeMinority = TRUE)
  mgp <-mgp_estimations$estimates %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
  rownames(count_df) <- count_df$Gene
  count_df <- count_df[,-1]
  bret_estimations<- BRETIGEA::findCells(count_df, bret_cell_markers, nMarker = 1000, method = "SVD",
                               scale = TRUE)
  bret <- bret_estimations %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Sample')
  mgp_final <- mgp
  bret_final <- bret
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
  compare_estimates <- merge(mgp, bret, by="Sample")
  compare_estimates <- compare_estimates[,-1]
  marker_mat <- data.matrix(compare_estimates, rownames.force = NA)
  cormat <- round(stats::cor(marker_mat),2)
  melted_cormat <- reshape::melt(cormat)
  names(melted_cormat)[3] <- "correlation"
  corr_plot <- ggplot2::ggplot(data = melted_cormat, ggplot2::aes(x=Var1, y=Var2, fill=correlation)) +
    ggplot2::geom_tile() + ggplot2::theme(axis.title.y=element_blank(),
                                          axis.title.x=element_blank(),
                                          axis.text.x = element_text(angle = 45, hjust = 1)) + ggplot2::geom_text(data=melted_cormat,aes(x=Var1, y=Var2,label=correlation)) +  ggplot2::scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar")

  return(list("mgp" = mgp_final, "bret" = bret_final, "corr_plot"= corr_plot))
}
