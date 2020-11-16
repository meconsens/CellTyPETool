#' Marker Genes for the Marker Gene Profile method
#'
#' Marker Genes from the BRETIGEA package reformatted for use in the Marker Gene Profile method
#'
#' @source BRETIGEA
#'
#' @format A nested list with 6 sub-lists of marker genes for brain cells: astrocyte, endothelial, microglia, neuron, oligodnedrocyte and oligodendrocyte precursor cells:
#' \describe{
#'  \item{ast}{Top 50 astrocyte marker genes.}
#'  \item{end}{Top 50 astrocyte marker genes.}
#'  \item{mic}{Top 50 microglia marker genes.}
#'  \item{neu}{Top 50 neuron marker genes.}
#'  \item{oli}{Top 50 oligodendrocyte marker genes.}
#'  \item{opc}{Top 50 oligodendrocyte precursor cell marker genes.}
#' }
#' @examples
#' \dontrun{
#'  mgp_cell_markers
#' }
"mgp_cell_markers"

#' Marker Genes for the BRETIGEA method
#'
#' Marker Genes from the BRETIGEA package
#'
#' @source BRETIGEA
#'
#' @format A dataframe with two columns: markers and cell.
#' \describe{
#'  \item{markers}{Top 1000 marker genes for  astrocyte, endothelial, microglia, neuron, oligodnedrocyte and oligodendrocyte precursor cells.}
#'  \item{cell}{Indicates which of the cell types ( astrocyte, endothelial, microglia, neuron, oligodnedrocyte and oligodendrocyte precursor cells) the marker gene is a marker for.}
#' }
#' @examples
#' \dontrun{
#'  bret_cell_markers
#' }
"bret_cell_markers"

#' Bulk-tissue RNAseq data counts
#'
#' Taken from the BRETIGEA package (originally aba_marker_expression)
#'
#' @source BRETIGEA
#'
#' @format A dataframe with 395 rows of genes and 378 subjects.
#' \describe{
#'  \item{Gene}{List of genes, some of which may be marker genes}
#'  \item{Subjects}{378 subjects}
#' }
#' @examples
#' \dontrun{
#'  count_df
#' }
"count_df"

#' Metadata for subjects corresponding to count_df subjects
#'
#' Taken from the BRETIGEA package, modified to have "Covariate" and "DiseasePhenotypeScore" columns from previous variables (names are not meaningful, just used to show relationships between cell type proportions and potentially a disease phenotype state, covarying for some confounding effect)
#'
#' @source BRETIGEA
#'
#' @format A dataframe with 345 subjects also in count_df and 3 rows indicating the subjects id (Sample), some covariate (covariate) and disease state score (DiseasePhenotypeScore)
#' \describe{
#'  \item{Sample}{The subjects id which corresponds to count_df}
#'  \item{Covariate}{Some confounding factor that must be accounted for when finding signals between cell-type proportions and a disease phenotype (could be batch etc.)}
#'  \item{DiseasePhenotypeScore}{Some indication of disease state, could be a score of combined AD pathology variables (CERAD score, Braak pathology stage etc.)}
#' }
#' @examples
#' \dontrun{
#'  metadata
#' }
"metadata"

# [END]
