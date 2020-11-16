
test_that("Returns list object", {

  data("countDf")
  data("bretCellMarkers")

  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]
  findCellsMod <- findCellsMod(
    inputMat = countDf,
    markers = bretCellMarkers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE)

  expect_type(findCellsMod, "list")
})

test_that("Returns items nested in list correctly",{
  data("countDf")
  data("bretCellMarkers")

  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]
  findCellsMod <- findCellsMod(
    inputMat = countDf,
    markers = bretCellMarkers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE)

  expect_type(findCellsMod$markerList, "list")
  expect_type(findCellsMod$SPVS, "double")
  expect_s3_class(findCellsMod$markerList, "data.frame")
})

context("Checking for invalid user input for findCellsMod")
test_that("findCellsMod error upon invalid user input", {

  data("countDf")
  data("bretCellMarkers")

  rownames(countDf) <- countDf$Gene
  countDf <- countDf[,-1]


  #no marker gene symbol present in the rownames of the input matrix
  wrongBretCellMarkers <- bretCellMarkers[,-1]
  expect_error(findCellsMod <- findCellsMod(
    inputMat = countDf,
    markers = wrongBretCellMarkers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE))

  # wrong method specified
  expect_error(findCellsMod <- findCellsMod(
    inputMat = countDf,
    markers = bretCellMarkers,
    nMarker = 10,
    method = "X",
    scale = TRUE))
})
# [END]

