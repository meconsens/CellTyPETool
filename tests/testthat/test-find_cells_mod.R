
test_that("Returns list object", {

  data("count_df")
  data("bret_cell_markers")

  rownames(count_df) <- count_df$Gene
  count_df <- count_df[,-1]
  find_cells_mod <- find_cells_mod(
    inputMat = count_df,
    markers = bret_cell_markers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE)

  expect_type(find_cells_mod, "list")
})

test_that("Returns items nested in list correctly",{
  data("count_df")
  data("bret_cell_markers")

  rownames(count_df) <- count_df$Gene
  count_df <- count_df[,-1]
  find_cells_mod <- find_cells_mod(
    inputMat = count_df,
    markers = bret_cell_markers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE)

  expect_type(find_cells_mod$marker_list, "list")
  expect_type(find_cells_mod$SPVS, "double")
  expect_s3_class(find_cells_mod$marker_list, "data.frame")
})

context("Checking for invalid user input for bretigea_marker_addition")
test_that("find_cells_mod error upon invalid user input", {

  data("count_df")
  data("bret_cell_markers")

  rownames(count_df) <- count_df$Gene
  count_df <- count_df[,-1]


  #no marker gene symbol present in the rownames of the input matrix
  wrong_bret_cell_markers <- bret_cell_markers[,-1]
  expect_error(find_cells_mod <- find_cells_mod(
    inputMat = count_df,
    markers = wrong_bret_cell_markers,
    nMarker = 10,
    method = "SVD",
    scale = TRUE))

  # wrong method specified
  expect_error(find_cells_mod <- find_cells_mod(
    inputMat = count_df,
    markers = wrong_bret_cell_markers,
    nMarker = 10,
    method = "X",
    scale = TRUE))
})
# [END]

