test_that("plot_spatial_gene returns a ggplot for different gene selectors", {
  spe <- mock_spe()

  p1 <- plot_spatial_gene(
    input = spe,
    gene = 1,
    assay_name = "counts",
    transform = "sqrt",
    point_size = 2
  )
  p2 <- plot_spatial_gene(
    input = spe,
    gene = "GENE3",
    assay_name = "counts",
    transform = "none"
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_spatial_gene validates gene and point size inputs", {
  spe <- mock_spe()

  expect_error(
    plot_spatial_gene(spe, gene = 99),
    "out of bounds"
  )
  expect_error(
    plot_spatial_gene(spe, gene = "missing_gene"),
    "was not found"
  )
  expect_error(
    plot_spatial_gene(spe, gene = 1, point_size = -1),
    "greater than or equal to 0"
  )
})

test_that("plot_spatial_fit returns a ggplot and supports covariates", {
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

  p <- plot_spatial_fit(
    input = spe,
    gene = "GENE3",
    assay_name = "counts",
    covariates = ~ batch,
    family = "poisson",
    use_bam = FALSE,
    smooth_k = 4,
    point_size = 2
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_spatial_fit errors when a single-gene fit cannot be obtained", {
  counts <- matrix(0, nrow = 4, ncol = 6)
  spe <- mock_spe(counts)

  expect_error(
    plot_spatial_fit(
      input = spe,
      gene = 1,
      assay_name = "counts",
      family = "poisson",
      use_bam = FALSE,
      smooth_k = 4
    ),
    "Library sizes could not be computed"
  )
})
