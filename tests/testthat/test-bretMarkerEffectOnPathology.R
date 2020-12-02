
test_that("Returns list object of ggplot type", {

  data("countDf")
  data("metadata")
  data("bretCellMarkers")

  bretMarkerEffectOnPathology <- bretMarkerEffectOnPathology(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore",
    cellTypeNames = unique(bretCellMarkers$cell),
    n= 10)

  expect_type(bretMarkerEffectOnPathology, "list")
  expect_s3_class(bretMarkerEffectOnPathology, c("gg", "ggplot"))
})

test_that("Scale labels are correct",{
  data("countDf")
  data("metadata")
  data("bretCellMarkers")

  bretMarkerEffectOnPathology <- bretMarkerEffectOnPathology(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathologyName = "DiseasePhenotypeScore",
    cellTypeNames = unique(bretCellMarkers$cell),
    n= 10)
  expect_identical(bretMarkerEffectOnPathology$labels$y, "-log10(p)")
  expect_identical(bretMarkerEffectOnPathology$labels$x, "nMarker")
})

context("Checking for invalid user input for bretMarkerEffectOnPathology")
test_that("bretMarkerEffectOnPathology error upon invalid user input", {

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
  expect_error(bretMarkerEffectOnPathology <- bretMarkerEffectOnPathology(
    countDf = countDf,
    bretCellMarkers = bretCellMarkers,
    cellName = "X",
    metadata = metadata,
    covar = covar,
    pathologyName = pathologyName,
    cellTypeNames = cellTypeNames,
    n= 10))

  # covar not in metadata
  expect_error(bretMarkerEffectOnPathology <- bretMarkerEffectOnPathology(
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
