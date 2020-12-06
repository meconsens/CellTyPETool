#' Launch Shiny App for CellTyPETool
#'
#' A function that launches the Shiny app for CellTyPETool
#' The purpose of this app is
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' CellTyPETool::runCellTyPETool()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runCellTyPETool <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "CellTyPETool")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
# [END]
