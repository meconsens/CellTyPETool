
test_that("Returns dataframe object", {

  data("countDf")
  data("mgpCellMarkers")

  mgpQCMetricsResults <- mgpQCMetrics(countDf, mgpCellMarkers)

  expect_s3_class(mgpQCMetricsResults, "data.frame")
})


context("Checking for invalid user input for mgpQCMetrics")
test_that("mgpQCMetrics error upon invalid user input", {

  data("countDf")
  data("mgpCellMarkers")

  #countDf doesn't have Gene column
  wrongCountDf <- countDf[,-1]

  # wrong method specified
  expect_error(mgpQCMetricsResults <- mgpQCMetrics(wrongCountDf, mgpCellMarkers))
})
# [END]

