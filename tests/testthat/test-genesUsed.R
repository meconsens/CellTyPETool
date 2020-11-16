
test_that("Returns list object", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  genesUsed <- genesUsed(
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  expect_type(genesUsed, "list")
})

test_that("Returns items nested in list correctly",{
  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  genesUsed <- genesUsed(
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  expect_type(genesUsed$mgpMarkers, "list")
  expect_type(genesUsed$bretMarkers, "list")
  expect_s3_class(genesUsed$mgpMarkers, "data.frame")
  expect_s3_class(genesUsed$bretMarkers, "data.frame")
})

context("Checking for invalid user input for genesUsed")
test_that("genesUsed error upon invalid user input", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  #no marker gene symbol present in the rownames of the input matrix
  wrongBretCellMarkers <- bretCellMarkers[,-1]
  expect_error(genesUsed <- genesUsed(
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = wrongBretCellMarkers))


  # countDf doesn't have Gene column
  wrongCountDf <- countDf[,-1]
  expect_error(genesUsed <- genesUsed(
    countDf = wrongCountDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers))
})
# [END]

