#' Calculates and compares the estimations of cell types derived by the
#' markerGeneProfile and BRETIGEA findCells method
#'
#' A function that calculates cell type proportions using two validated methods,
#' markerGeneProfile and findCells, and then compares the cell type proportions
#' calculated by each
#'
#' @param estimates The estimates of cell type proportions returned by calc_and_compare()
#' either markerGeneProfile derived or BRETIGEA derived
#' @param metadata A dataframe with subjects also in count_df and rows
#' indicating the subjects id, some covariate, and disease state score or
#' pathology.
#' @param cell_type_names The names of all the unique cell types for which
#' there are markers in bret_cell_markers: unique(bret_cell_markers$cell)
#' @param covar A covariate to be taken into account when running linear models
#' to check the association between the cell type indicated by cell and the
#' pathology indicated by pathology_name.
#' @param pathology_name The pathology associated with the disease in question
#' for which the association between it and the cell type indicated by cell is
#' being examined.
#'
#' @return Returns a volcano plot showing the significance of assocations between
#' each cell type proportion derived and the pathology in question
#'
#' @examples
#' # Examples 1:
#' # Using count_df, bret_cell_markers, mgp_cell_markers datasets available with package
#'
#'calc_and_compare <- calc_and_compare (
#'                 count_df = count_df,
#'                 mgp_cell_markers = mgp_cell_markers,
#'                 bret_cell_markers = bret_cell_markers)
#'
#' estimates_v_phenotype <- estimates_v_phenotype(
#'                 estimates = calc_and_compare$bret,
#'                 metadata = metadata,
#'                 cell_type_names = unique(bret_cell_markers$cell),
#'                 covar = "Covariate",
#'                 pathology = "DiseasePhenotypeScore")
#'
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
#' York, 2016.
#'
#' Kamil Slowikowski (2020). ggrepel: Automatically Position Non-Overlapping Text
#' Labels with 'ggplot2'. R package version 0.8.2.
#' https://CRAN.R-project.org/package=ggrepel
#'
#' Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., Rocco, B., Sibille, E., & Pavlidis, P. (2017).
#' CrossLaboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk
#' Tissue Data. \emph{eNeuro}, 4(6), ENEURO.0212-17.2017. https://doi.org/10.1523/ENEURO.0212-17.201
#'
#' McKenzie, A.T., Wang, M., Hauberg, M.E. et al. Brain Cell Type Specific Gene Expression and Coexpression Network Architectures.
#' \emph{Sci Rep} 8, 8868 (2018). https://doi.org/10.1038/s41598-018-27293-5
#'
#'R Core Team (2020). R: A language and environment for statistical computing. R
#'Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#'
#' Wickham et al., (2019). Welcome to the tidyverse.\emph{ Journal of Open Source Software}, 4(43), 1686,
#' https://doi.org/10.21105/joss.01686
#'
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#' @export
#' @importFrom stats lm anova coef p.adjust as.formula
#' @import ggplot2
#' @import ggrepel
estimates_v_phenotype <-function(estimates, metadata, cell_type_names, covar, pathology_name){

  if(!covar %in% colnames(metadata)) stop("The covar argument must have a corresponding column in metadata.")

  if(!pathology_name %in% colnames(metadata)) stop("The pathology_name argument must have a corresponding column in metadata.")


  model.data <- merge(estimates, metadata)
  results <- sapply(cell_type_names,function(celltype) {
    sapply(pathology_name, function(pathology) {
      form <- as.formula(paste0(pathology,"~",celltype," + ",paste0(covar,collapse=" + ")))
      model <- stats::lm(data=model.data,form)
      p <- stats::anova(model)[celltype,5]
      beta <- stats::coef(model)[celltype]
      n <- nrow(model$model)
      return(c(celltype,pathology,beta, p,n))

    })
  })

  results <- as.data.frame(matrix(results,ncol=5,byrow = T),stringsAsFactors = F)
  names(results) <- c("celltype","pathology","beta","p","n")

  results <- within(results,{
    p <- as.numeric(p)
    beta <- as.numeric(beta)
    n <- as.numeric(n)})

  results$bonfp <- stats::p.adjust(results$p, method="bonferroni")
  results$fdr <- stats:: p.adjust(results$p, method="fdr")

  volcano_plot <- ggplot2::ggplot(data=results,ggplot2::aes(y=-log10(fdr),x=beta,col=celltype))+
    ggplot2::geom_hline(yintercept = -log10(0.05),col="blue",lty=3)+
    ggrepel::geom_text_repel(data=subset(results,p<0.05),ggplot2::aes(label=celltype))+
    ggplot2::geom_vline(xintercept = 0,col="red")+ggplot2::geom_point(size=1.5)+
    ggplot2::facet_wrap(~pathology)+ ggplot2::theme_minimal() +ggplot2::theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
  return(volcano_plot)

}
