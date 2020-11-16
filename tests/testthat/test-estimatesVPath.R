
test_that("Returns list object of ggplot type", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")


  calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  estimatesVPath <- estimatesVPath(
    estimates = calcAndCompare$bret,
    metadata = metadata,
    cellTypeNames = unique(bretCellMarkers$cell),
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore")


  expect_type(estimatesVPath, "list")
  expect_s3_class(estimatesVPath, c("gg", "ggplot"))
})

test_that("Scale labels are correct",{
  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")


  calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  estimatesVPath <- estimatesVPath(
    estimates = calcAndCompare$bret,
    metadata = metadata,
    cellTypeNames = unique(bretCellMarkers$cell),
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore")

  expect_identical(estimatesVPath$labels$y, "-log10(fdr)")
  expect_identical(estimatesVPath$labels$x, "beta")
})

context("Checking for invalid user input for estimatesVPath")
test_that("estimatesVPath error upon invalid user input", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  covar <- "Covariate"
  pathologyName <- "DiseasePhenotypeScore"


  calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)



  # covar not in metadata
  expect_error(estimatesVPath <- estimatesVPath(
    estimates = calcAndCompare$bret,
    metadata = metadata,
    cellTypeNames = unique(bretCellMarkers$cell),
    covar = "X",
    pathologyName = "DiseasePhenotypeScore"))

  # pathology not in metadata
  expect_error(  estimatesVPath <- estimatesVPath(
    estimates = calcAndCompare$bret,
    metadata = metadata,
    cellTypeNames = unique(bretCellMarkers$cell),
    covar = "Covariate",
    pathologyName = "X"))
})
# [END]
