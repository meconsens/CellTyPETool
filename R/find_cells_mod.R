#' Modified BRETIGEA findCells() function
#'
#' A function that runs BRETIGEA findCells() but returns the used markers
#'
#' @param inputMat Numeric gene expression data frame or matrix, with
#' rownames corresponding to gene names, some of which are marker genes,
#' and columns corresponding to samples.
#' @param markers Data frame with marker genes in one column (named "marker")
#' and the cell type that that gene symbol corresponds to in another column
#' (named "cell").
#' @param nMarker The number of marker genes (that are present in your expression
#'  data set) to use in estimating the surrogate cell type proportion variable
#'  for each cell type.
#' @param method 	To estimate the cell type proportions, can either use PCA or SVD.
#' @param scale Whether or not to scale the gene expression data from each
#' marker gene prior to using it as an input for dimension reduction.
#'
#' @return Returns a  list of sample-by-cell type matrix of estimate cell
#' type proportionvariables and the markers used
#'
#' @examples
#' # Examples 1:
#' # Using bret_cell_markers, metadata datasets available with package
#'
#' rownames(count_df) <- count_df$Gene
#' count_df <- count_df[,-1]
#' find_cells_mod <- find_cells_mod(
#'                inputMat = count_df,
#'                markers = bret_cell_markers,
#'                nMarker = 10,
#'                method = "SVD",
#'                scale = TRUE)
#' find_cells_mod$marker_list
#' find_cells_mod$SPVS
#'
#' @references
#' #' Chikina M, Zaslavsky E, Sealfon SC. CellCODE: a robust latent variable approach to differential expression analysis for heterogeneous cell populations. \emph{Bioinformatics}
#' . 2015;31(10):1584-91.
#'
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
#' URL https://www.R-project.org/.
#' @export
#' @importFrom stats cor prcomp
find_cells_mod <- function(inputMat, markers, nMarker, method, scale){
  if(!all(c("markers", "cell") %in% colnames(markers))){
    stop("The markers argument must be a df with a column named marker s(gene symbols) and a column named cell (corresponding cell types).")
  }

  if(!any(markers$markers %in% rownames(inputMat))) stop("At least one marker gene symbol must be present in the rownames of the input matrix.")

  if(!method %in% c("SVD", "PCA")) stop("The method argument must be either SVD or PCA.")

  cell_types = unique(markers$cell)

  #choose the appropriate number of top markers from each cell type that are actually present in the data set
  for(i in 1:length(cell_types)){
    top_markers_tmp = markers[markers$cell == cell_types[i], ]
    top_markers_present = vector()
    top_markers_present_count = 0
    for(j in 1:nrow(top_markers_tmp)){
      gene = top_markers_tmp[j, "markers"]
      if(gene %in% rownames(inputMat)){
        if(sum(inputMat[gene, ] > 0)){
          top_markers_present = c(top_markers_present, gene)
          top_markers_present_count = top_markers_present_count + 1
        }
      }
      if(top_markers_present_count == nMarker) break
    }
    tmp_markers_top = data.frame(markers = top_markers_present,
                                 cell = rep(cell_types[i], length(top_markers_present)), stringsAsFactors = FALSE)
    if(i == 1){
      markers_df = tmp_markers_top
    }
    if(i > 1){
      markers_df = rbind(markers_df, tmp_markers_top)
    }
  }


  #adapted from CellCODE
  for(i in 1:length(cell_types)){
    genes = markers_df[markers_df$cell == cell_types[i], "markers"]
    data = inputMat[genes, ]
    if(scale){
      means = apply(data, 1, mean)
      data = sweep(data, 1, means, "-")
    }
    if(method == "PCA"){
      pcres = stats::prcomp(t(data))
      props = pcres$x[ , 1]
    }
    if(method == "SVD"){
      svdres = svd(data)
      props = svdres$v[ , 1]
    }
    cor_res = stats::cor(props, t(data[genes, ]))
    if (mean(cor_res, na.rm = TRUE) < 0) {
      props = -props
    }
    if(i == 1){
      SPVs = props
    } else {
      SPVs = cbind(SPVs, props)
    }
  }

  colnames(SPVs) = cell_types
  rownames(SPVs) = colnames(inputMat)

  markersandcells <- list("marker_list" = markers_df, "SPVS" = SPVs)
  return(markersandcells)
}
