#' Displays strength of association between pathologies and cell types
#'
#' A function that generates a volcano plot showing the significance of associations between
#' each cell type proportion derived and the pathology in question
#'
#' @param estimates The estimates of cell type proportions returned by calcAndCompare()
#'     either markerGeneProfile derived or BRETIGEA derived
#' @param metadata A dataframe with subjects also in countDf and rows
#'     indicating the subjects id, some covariate, and disease state score or
#'     pathology.
#' @param cellTypeNames The names of all the unique cell types for which
#'     there are markers in bretCellMarkers: unique(bretCellMarkers$cell)
#' @param covar A covariate to be taken into account when running linear models
#'     to check the association between the cell type indicated by cell and the
#'     pathology indicated by pathologyName.
#' @param pathologyName The pathology associated with the disease in question
#'     for which the association between it and the cell type indicated by cell is
#'     being examined.
#'
#' @return Returns a volcano plot showing the significance of associations between
#' each cell type proportion derived and the pathology in question
#'
#' @examples
#' # Examples 1:
#' # Using countDf, bretCellMarkers, mgpCellMarkers datasets available with package
#'
#'calcAndCompareResults <- calcAndCompare (
#'                 countDf = countDf,
#'                 mgpCellMarkers = mgpCellMarkers,
#'                 bretCellMarkers = bretCellMarkers)
#'
#' estimatesVPathResults <- estimatesVPath(
#'                 estimates = calcAndCompareResults$bret,
#'                 metadata = metadata,
#'                 cellTypeNames = unique(bretCellMarkers$cell),
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
estimatesVPath <-function(estimates, metadata, cellTypeNames, covar, pathologyName){

  if(!covar %in% colnames(metadata)) stop("The covar argument must have a corresponding column in metadata.")

  if(!pathologyName %in% colnames(metadata)) stop("The pathologyName argument must have a corresponding column in metadata.")

  #calculate direction of effect and significance of association between cell type and pathology, as well as number of observations
  model.data <- merge(estimates, metadata)
  results <- sapply(cellTypeNames,function(celltype) {
    sapply(pathologyName, function(pathology) {
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
  #correctly assign types (i.e. p value is a numeric value)
  results <- within(results,{
    p <- as.numeric(p)
    beta <- as.numeric(beta)
    n <- as.numeric(n)})

  #correct p value fior multiple tests
  results$bonfp <- stats::p.adjust(results$p, method="bonferroni")
  results$fdr <- stats:: p.adjust(results$p, method="fdr")

  #generate volcano plot
  volcanoPlot <- ggplot2::ggplot(data=results,ggplot2::aes(y=-log10(fdr),x=beta,col=celltype))+
    ggplot2::geom_hline(yintercept = -log10(0.05),col="blue",lty=3)+
    ggrepel::geom_text_repel(data=subset(results,p<0.05),ggplot2::aes(label=celltype))+
    ggplot2::geom_vline(xintercept = 0,col="red")+ggplot2::geom_point(size=1.5)+
    ggplot2::facet_wrap(~pathology)+ ggplot2::theme_minimal() +ggplot2::theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
  return(volcanoPlot)

}
#[END]
