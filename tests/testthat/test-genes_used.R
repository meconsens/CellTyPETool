
test_that("Returns list object", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  genes_used <- genes_used(
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

  expect_type(genes_used, "list")
})

test_that("Returns items nested in list correctly",{
  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  genes_used <- genes_used(
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

  expect_type(genes_used$mgp_markers, "list")
  expect_type(genes_used$bret_markers, "list")
  expect_s3_class(genes_used$mgp_markers, "data.frame")
  expect_s3_class(genes_used$bret_markers, "data.frame")
})

context("Checking for invalid user input for bretigea_marker_addition")
test_that("genes_used error upon invalid user input", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  #no marker gene symbol present in the rownames of the input matrix
  wrong_bret_cell_markers <- bret_cell_markers[,-1]
  expect_error(genes_used <- genes_used(
      count_df = count_df,
      mgp_cell_markers = mgp_cell_markers,
      bret_cell_markers = wrong_bret_cell_markers))


  # count_df doesn't have Gene column
  wrong_count_df <- count_df[,-1]
  expect_error(genes_used <- genes_used(
      count_df = wrong_count_df,
      mgp_cell_markers = mgp_cell_markers,
      bret_cell_markers = bret_cell_markers))
})
# [END]

