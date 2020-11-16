
test_that("Returns list object of ggplot type", {

  data("count_df")
  data("metadata")
  data("bret_cell_markers")

  bretigea_marker_addition <- bretigea_marker_addition(
    count_df = count_df,
    bret_cell_markers = bret_cell_markers,
    cell_name = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathology_name = "DiseasePhenotypeScore",
    cell_type_names = unique(bret_cell_markers$cell),
    n= 10)

  expect_type(bretigea_marker_addition, "list")
  expect_s3_class(bretigea_marker_addition, c("gg", "ggplot"))
})

test_that("Scale labels are correct",{
  data("count_df")
  data("metadata")
  data("bret_cell_markers")

  bretigea_marker_addition <- bretigea_marker_addition(
    count_df = count_df,
    bret_cell_markers = bret_cell_markers,
    cell_name = "mic",
    metadata = metadata,
    covar = "Covariate",
    pathology_name = "DiseasePhenotypeScore",
    cell_type_names = unique(bret_cell_markers$cell),
    n= 10)
  expect_identical(bretigea_marker_addition$labels$y, "-log10(p)")
  expect_identical(bretigea_marker_addition$labels$x, "nMarker")
})

context("Checking for invalid user input for bretigea_marker_addition")
test_that("bretigea_marker_addition error upon invalid user input", {

  data("count_df")
  data("metadata")
  data("bret_cell_markers")
  count_df = count_df
  bret_cell_markers = bret_cell_markers
  cell_name = "mic"
  metadata = metadata
  covar = "Covariate"
  pathology_name = "DiseasePhenotypeScore"
  cell_type_names = unique(bret_cell_markers$cell)
  n= 10

  # cell_name not in cell_type_names
  expect_error(bretigea_marker_addition <- bretigea_marker_addition(
    count_df = count_df,
    bret_cell_markers = bret_cell_markers,
    cell_name = "X",
    metadata = metadata,
    covar = covar,
    pathology_name = pathology_name,
    cell_type_names = cell_type_names,
    n= 10))

  # covar not in metadata
  expect_error(bretigea_marker_addition <- bretigea_marker_addition(
    count_df = count_df,
    bret_cell_markers = bret_cell_markers,
    cell_name = cell_name,
    metadata = metadata,
    covar = "X",
    pathology_name = pathology_name,
    cell_type_names = cell_type_names,
    n= 10))
})
# [END]
