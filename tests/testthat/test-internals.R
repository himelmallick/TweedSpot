test_that("tweedspot_libsize uses sizeFactor when present and repairs bad values", {
  spe <- mock_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  SummarizedExperiment::colData(spe)$sizeFactor <- c(1, NA, 2, 0, 3, 4)

  libsz <- TweedSpot:::tweedspot_libsize(spe, Y)

  expect_equal(libsz, c(1, 1, 2, 1, 3, 4))
})

test_that("tweedspot_check_scalar validates numeric bounds", {
  expect_no_error(TweedSpot:::tweedspot_check_scalar(1, "x", min_value = 0, max_value = 2))
  expect_error(TweedSpot:::tweedspot_check_scalar("a", "x"), "single numeric value")
  expect_error(TweedSpot:::tweedspot_check_scalar(-1, "x", min_value = 0), "greater than or equal to 0")
  expect_error(TweedSpot:::tweedspot_check_scalar(3, "x", max_value = 2), "less than or equal to 2")
})

test_that("tweedspot_order_results sorts with statistic tie-breaking", {
  res <- data.frame(
    tweedspot_padj = c(0.1, 0.1, NA),
    tweedspot_stat = c(2, 5, 1)
  )

  asc <- TweedSpot:::tweedspot_order_results(res, key_col = "tweedspot_padj", decreasing = FALSE)
  desc <- TweedSpot:::tweedspot_order_results(res, key_col = "tweedspot_stat", decreasing = TRUE)

  expect_equal(asc$tweedspot_stat, c(5, 2, 1))
  expect_equal(desc$tweedspot_stat, c(5, 2, 1))
})

test_that("annotation exclusions support name and biotype based filtering", {
  spe <- mock_spe()
  rd <- SummarizedExperiment::rowData(spe)
  rd$biotype <- c("mitochondrial_gene", "ribosomal protein", "other", "other")
  SummarizedExperiment::rowData(spe) <- rd

  excluded <- TweedSpot:::tweedspot_annotation_exclusions(
    input = spe,
    exclude_mito = TRUE,
    exclude_ribo = TRUE,
    gene_name_col = "gene_name",
    gene_biotype_col = "biotype"
  )

  expect_identical(excluded, c(TRUE, TRUE, FALSE, FALSE))
  expect_error(
    TweedSpot:::tweedspot_annotation_exclusions(
      input = spe,
      exclude_mito = TRUE,
      exclude_ribo = FALSE,
      gene_name_col = "bad_col"
    ),
    "gene_name_col"
  )
  expect_error(
    TweedSpot:::tweedspot_annotation_exclusions(
      input = spe,
      exclude_mito = FALSE,
      exclude_ribo = TRUE,
      gene_biotype_col = "bad_col"
    ),
    "gene_biotype_col"
  )
})

test_that("covariate parsing supports formulas and character vectors", {
  spe <- mock_spe()
  SummarizedExperiment::colData(spe)$group_chr <- rep(c("x", "y"), length.out = ncol(spe))

  cov_formula <- TweedSpot:::tweedspot_covariates(spe, ~ batch)
  cov_char <- TweedSpot:::tweedspot_covariates(spe, "group_chr")

  expect_true(is.list(cov_formula))
  expect_true(is.factor(cov_char$data$group_chr))
  expect_match(cov_char$terms, "group_chr")
})

test_that("covariate parsing rejects invalid inputs", {
  spe <- mock_spe()

  expect_null(TweedSpot:::tweedspot_covariates(spe, NULL))
  expect_error(TweedSpot:::tweedspot_covariates(spe, ~ 1), "must reference at least one column")
  expect_error(TweedSpot:::tweedspot_covariates(spe, "missing"), "not found")

  SummarizedExperiment::colData(spe)$bad <- I(vector("list", ncol(spe)))
  expect_error(TweedSpot:::tweedspot_covariates(spe, "bad"), "must be numeric, logical, character, or factor")

  spe_na <- mock_spe()
  SummarizedExperiment::colData(spe_na)$batch[1] <- NA
  expect_error(TweedSpot:::tweedspot_covariates(spe_na, ~ batch), "must not contain missing values")

  spe_reserved <- mock_spe()
  SummarizedExperiment::colData(spe_reserved)$expr <- seq_len(ncol(spe_reserved))
  expect_error(TweedSpot:::tweedspot_covariates(spe_reserved, ~ expr), "reserved names")
})

test_that("gene resolution and coordinate helpers validate inputs", {
  spe <- mock_spe()

  resolved_by_index <- TweedSpot:::tweedspot_resolve_gene(spe, 1)
  resolved_by_name <- TweedSpot:::tweedspot_resolve_gene(spe, "GENE3")
  coords <- TweedSpot:::tweedspot_plot_coords(spe)

  expect_equal(resolved_by_index$index, 1L)
  expect_match(resolved_by_name$label, "GENE3")
  expect_equal(ncol(coords), 2)

  expect_error(TweedSpot:::tweedspot_resolve_gene(spe, c(1, 2)), "must be a single row index")
})

test_that("basis selection and cct helper cover their branches", {
  coords_small <- cbind(x = c(1, 1, 2), y = c(1, 1, 2))
  coords_large <- cbind(x = 1:6, y = c(1, 2, 3, 1, 2, 3))

  expect_equal(TweedSpot:::tweedspot_basis_k(coords_small), 3L)
  expect_equal(TweedSpot:::tweedspot_basis_k(coords_large, smooth_k = 100), 5L)
  expect_equal(TweedSpot:::tweedspot_basis_k(coords_large, smooth_k = 2), 3L)

  expect_true(is.finite(TweedSpot:::cct(c(0.2, 0.4))))
  expect_true(is.finite(TweedSpot:::cct(c(1e-20, 0.4), weights = c(2, 1))))
})

test_that("single-gene fitting helper and agnostic engine cover major branches", {
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
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  coords <- scale(SpatialExperiment::spatialCoords(spe))
  libsz <- TweedSpot:::tweedspot_libsize(spe, Y)
  covariates <- TweedSpot:::tweedspot_covariates(spe, ~ batch)

  fit <- TweedSpot:::tweedspot_fit_single_gene(
    expr = as.numeric(Y[1, ]),
    coords = coords,
    libsz = libsz,
    covariates = covariates,
    family = "poisson",
    fit_method = "REML",
    use_bam = FALSE,
    bam_discrete = TRUE,
    bam_nthreads = 1L,
    smooth_k = 4
  )

  expect_s3_class(fit, "gam")

  one_part <- suppressWarnings(
    TweedSpot:::tweedspot_agnostic(
      Y = Y,
      coords = coords,
      libsz = libsz,
      covariates = covariates,
      two_part = FALSE,
      combine = "CCT",
      family = "poisson",
      fit_method = "REML",
      use_bam = FALSE,
      bam_discrete = TRUE,
      bam_nthreads = 1L,
      smooth_k = 4,
      BPPARAM = BiocParallel::SerialParam()
    )
  )

  two_part <- suppressWarnings(
    TweedSpot:::tweedspot_agnostic(
      Y = Y,
      coords = coords,
      libsz = libsz,
      covariates = covariates,
      two_part = TRUE,
      combine = "stouffer",
      family = "poisson",
      fit_method = "REML",
      use_bam = FALSE,
      bam_discrete = TRUE,
      bam_nthreads = 1L,
      smooth_k = 4,
      BPPARAM = BiocParallel::SerialParam()
    )
  )

  min_part <- suppressWarnings(
    TweedSpot:::tweedspot_agnostic(
      Y = Y,
      coords = coords,
      libsz = libsz,
      covariates = covariates,
      two_part = TRUE,
      combine = "min",
      family = "poisson",
      fit_method = "REML",
      use_bam = FALSE,
      bam_discrete = TRUE,
      bam_nthreads = 1L,
      smooth_k = 4,
      BPPARAM = BiocParallel::SerialParam()
    )
  )

  expect_equal(length(one_part$pval), nrow(Y))
  expect_equal(length(two_part$pval), nrow(Y))
  expect_equal(length(min_part$pval), nrow(Y))

  expect_error(
    TweedSpot:::tweedspot_fit_single_gene(
      expr = as.numeric(Y[1, ]),
      coords = coords,
      libsz = c(0, rep(1, ncol(Y) - 1)),
      covariates = covariates,
      family = "poisson",
      fit_method = "REML",
      use_bam = FALSE,
      bam_discrete = TRUE,
      bam_nthreads = 1L,
      smooth_k = 4
    ),
    "positive and finite"
  )

  expect_error(
    TweedSpot:::tweedspot_agnostic(
      Y = Y,
      coords = coords,
      libsz = c(0, rep(1, ncol(Y) - 1)),
      covariates = covariates,
      two_part = FALSE,
      combine = "CCT",
      family = "poisson",
      fit_method = "REML",
      use_bam = FALSE,
      bam_discrete = TRUE,
      bam_nthreads = 1L,
      smooth_k = 4,
      BPPARAM = BiocParallel::SerialParam()
    ),
    "positive and finite"
  )
})
