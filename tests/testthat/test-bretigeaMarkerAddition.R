
test_that("Returns list object of ggplot type", {

  data("countDf")
  data("metadata")
  data("bretCellMarkers")

  bretigeaMarkerAddition <- bretigeaMarkerAddition(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore",
    cellTypeNames = unique(bretCellMarkers$cell),
    n= 10)

  expect_type(bretigeaMarkerAddition, "list")
  expect_s3_class(bretigeaMarkerAddition, c("gg", "ggplot"))
})

test_that("Scale labels are correct",{
  data("countDf")
  data("metadata")
  data("bretCellMarkers")

  bretigeaMarkerAddition <- bretigeaMarkerAddition(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore",
    cellTypeNames = unique(bretCellMarkers$cell),
    n= 10)
  expect_identical(bretigeaMarkerAddition$labels$y, "-log10(p)")
  expect_identical(bretigeaMarkerAddition$labels$x, "nMarker")
})

context("Checking for invalid user input for bretigeaMarkerAddition")
test_that("bretigeaMarkerAddition error upon invalid user input", {

  data("countDf")
  data("metadata")
  data("bretCellMarkers")
  countDf = countDf
  bretCellMarkers = bretCellMarkers
  cellName = "mic"
  metadata = metadata
  covar = "Covariate"
  pathologyName = "DiseasePhenotypeScore"
  cellTypeNames = unique(bretCellMarkers$cell)
  n= 10

  # cellName not in cellTypeNames
  expect_error(bretigeaMarkerAddition <- bretigeaMarkerAddition(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "X",
    metadata = metadata,
    covar = covar,
    pathologyName = pathologyName,
    cellTypeNames = cellTypeNames,
    n= 10))

  # covar not in metadata
  expect_error(bretigeaMarkerAddition <- bretigeaMarkerAddition(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = cellName,
    metadata = metadata,
    covar = "X",
    pathologyName = pathologyName,
    cellTypeNames = cellTypeNames,
    n= 10))
})
# [END]
