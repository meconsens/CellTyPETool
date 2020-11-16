
test_that("Returns list object of ggplot type", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")


  calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

 estimates_v_phenotype <- estimates_v_phenotype(
   estimates = calc_and_compare$bret,
   metadata = metadata,
   cell_type_names = unique(bret_cell_markers$cell),
   covar = "Covariate",
   pathology_name = "DiseasePhenotypeScore")


 expect_type(estimates_v_phenotype, "list")
 expect_s3_class(estimates_v_phenotype, c("gg", "ggplot"))
})

test_that("Scale labels are correct",{
  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")


  calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)

  estimates_v_phenotype <- estimates_v_phenotype(
    estimates = calc_and_compare$bret,
    metadata = metadata,
    cell_type_names = unique(bret_cell_markers$cell),
    covar = "Covariate",
    pathology_name = "DiseasePhenotypeScore")

  expect_identical(estimates_v_phenotype$labels$y, "-log10(fdr)")
  expect_identical(estimates_v_phenotype$labels$x, "beta")
})

context("Checking for invalid user input for bretigea_marker_addition")
test_that("estimates_v_phenotype error upon invalid user input", {

  data("count_df")
  data("mgp_cell_markers")
  data("bret_cell_markers")

  covar <- "Covariate"
  pathology_name <- "DiseasePhenotypeScore"


  calc_and_compare <- calc_and_compare (
    count_df = count_df,
    mgp_cell_markers = mgp_cell_markers,
    bret_cell_markers = bret_cell_markers)



  # covar not in metadata
  expect_error(estimates_v_phenotype <- estimates_v_phenotype(
    estimates = calc_and_compare$bret,
    metadata = metadata,
    cell_type_names = unique(bret_cell_markers$cell),
    covar = "X",
    pathology_name = "DiseasePhenotypeScore"))

  # pathology not in metadata
  wrong_count_df <- count_df[,-1]
  expect_error(  estimates_v_phenotype <- estimates_v_phenotype(
    estimates = calc_and_compare$bret,
    metadata = metadata,
    cell_type_names = unique(bret_cell_markers$cell),
    covar = "Covariate",
    pathology_name = "X"))
})
# [END]
