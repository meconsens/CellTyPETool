
test_that("Returns list object", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

  expect_type(calc_and_compare, "list")
})

test_that("Returns items nested in list correctly",{
  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

  expect_type(calc_and_compare$mgp, "list")
  expect_type(calc_and_compare$bret, "list")
  expect_type(calc_and_compare$corr_plot, "list")
  expect_s3_class(calc_and_compare$mgp, "data.frame")
  expect_s3_class(calc_and_compare$bret, "data.frame")
})

context("Checking for invalid user input for bretigea_marker_addition")
test_that("calc_and_compare error upon invalid user input", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")



  # wrong bret_cell_markers format
  wrong_bret_cell_markers <- bret_cell_markers
  colnames(wrong_bret_cell_markers) <- c("wrong", "cell")
  expect_error(calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = wrong_bret_cell_markers))

  # count_df doesn't have Gene column
  wrong_count_df <- count_df[,-1]
  expect_error(calc_and_compare <- calc_and_compare (
    count_df = wrong_count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers))
})
# [END]
