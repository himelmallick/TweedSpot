test_that("tweedspot returns result columns for a small dataset", {
  counts <- matrix(
    c(
      0, 1, 2, 3, 4, 5, 5, 4,
      5, 4, 3, 2, 1, 0, 0, 1,
      1, 1, 1, 2, 2, 2, 3, 3,
      0, 0, 1, 1, 2, 2, 3, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  spe <- mock_spe(counts)

  fit <- tweedspot(
    input = spe,
    assay_name = "counts",
    covariates = ~ batch,
    two_part = FALSE,
    family = "poisson",
    use_bam = FALSE,
    smooth_k = 4,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )

  rd <- SummarizedExperiment::rowData(fit)
  expect_true(all(c(
    "tweedspot_stat",
    "tweedspot_pval",
    "tweedspot_padj",
    "tweedspot_edf",
    "tweedspot_dev_expl"
  ) %in% colnames(rd)))
})

test_that("tweedspot rejects constant covariates", {
  spe <- mock_spe()
  SummarizedExperiment::colData(spe)$constant <- "A"

  expect_error(
    tweedspot(
      input = spe,
      assay_name = "counts",
      covariates = ~ constant,
      family = "poisson",
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    ),
    "must vary across spatial locations"
  )
})

test_that("tweedspot errors on non-positive library sizes", {
  counts <- matrix(0, nrow = 3, ncol = 4)
  spe <- mock_spe(counts)

  expect_error(
    tweedspot(
      input = spe,
      assay_name = "counts",
      family = "poisson",
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    ),
    "Library sizes could not be computed"
  )
})
