
test_that("Returns list object", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  expect_type(calcAndCompare, "list")
})

test_that("Returns items nested in list correctly",{
  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")

  calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers)

  expect_type(calcAndCompare$mgp, "list")
  expect_type(calcAndCompare$bret, "list")
  expect_type(calcAndCompare$corrPlot, "list")
  expect_s3_class(calcAndCompare$mgp, "data.frame")
  expect_s3_class(calcAndCompare$bret, "data.frame")
})

context("Checking for invalid user input for calcAndCompare")
test_that("calcAndCompare error upon invalid user input", {

  data("countDf")
  data("mgpCellMarkers")
  data("bretCellMarkers")



  # wrong bretCellMarkers format
  wrongBretCellMarkers <- bretCellMarkers
  colnames(wrongBretCellMarkers) <- c("wrong", "cell")
  expect_error(calcAndCompare <- calcAndCompare (
    countDf = countDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = wrongBretCellMarkers))

  # countDf doesn't have Gene column
  wrongCountDf <- countDf[,-1]
  expect_error(calcAndCompare <- calcAndCompare (
    countDf = wrongCountDf,
    mgpCellMarkers = mgpCellMarkers,
    bretCellMarkers = bretCellMarkers))
})
# [END]
