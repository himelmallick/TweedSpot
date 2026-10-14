# Assays stored as sparse or DelayedArray matrices, and coordinate handling.

storage_spe <- function(assay) {
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = assay),
    spatialCoords = as.matrix(expand.grid(x = 1:12, y = 1:12))
  )
}

storage_counts <- function() {
  set.seed(8)
  grid <- expand.grid(x = 1:12, y = 1:12)
  counts <- rbind(
    band = stats::rpois(144, ifelse(grid$y > 6, 10, 1)),
    sparse = stats::rbinom(144, 3, 0.15),
    flat = stats::rpois(144, 4)
  )
  colnames(counts) <- paste0("loc", seq_len(144))
  counts
}

test_that("sparse and DelayedArray assays give the same results as dense ones", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("DelayedArray")
  counts <- storage_counts()
  run <- function(assay) {
    SummarizedExperiment::rowData(tweedspot(
      storage_spe(assay), smooth_k = 8, size_factor_method = "library_size",
      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
    ))$tweedspot_pval
  }
  dense <- run(counts)
  expect_equal(run(methods::as(counts, "CsparseMatrix")), dense, tolerance = 1e-8)
  expect_equal(run(DelayedArray::DelayedArray(counts)), dense, tolerance = 1e-8)

  flt <- function(assay) rownames(filter_features(storage_spe(assay), min_total = 0,
                                                   min_var = 0.5, smooth_k = 8))
  expect_identical(flt(methods::as(counts, "CsparseMatrix")), flt(counts))
  expect_identical(flt(DelayedArray::DelayedArray(counts)), flt(counts))
})

test_that("sparse-safe helpers agree with their dense versions", {
  skip_if_not_installed("Matrix")
  counts <- storage_counts()
  sp <- methods::as(counts, "CsparseMatrix")
  expect_equal(tweedspot_row_vars(sp), unname(apply(counts, 1, stats::var)))
  expect_equal(tweedspot_row_vars(counts), unname(apply(counts, 1, stats::var)))
  expect_true(tweedspot_is_integer_valued(sp))
  expect_false(tweedspot_is_integer_valued(sp * 0.5 + 0.25))
  skip_if_not_installed("DelayedArray")
  expect_true(tweedspot_is_integer_valued(DelayedArray::DelayedArray(counts)))
  expect_false(tweedspot_is_integer_valued(DelayedArray::DelayedArray(counts + 0.5)))
})

test_that("missing values in the assay are rejected clearly", {
  counts <- storage_counts()
  counts[1, 1] <- NA
  expect_error(filter_features(storage_spe(counts)), "missing values")
  expect_error(tweedspot(storage_spe(counts), verbose = FALSE), "missing values")
})

test_that("coordinates share one scale, so the tissue's shape is kept", {
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(1, 2, 60)),
    spatialCoords = cbind(x = rep(seq(0, 58, by = 2), 2), y = rep(c(0, 6), each = 30))
  )
  xy <- tweedspot_model_coords(spe)
  raw <- SpatialExperiment::spatialCoords(spe)
  expect_equal(diff(range(xy[, 1])) / diff(range(xy[, 2])),
               diff(range(raw[, 1])) / diff(range(raw[, 2])))
  expect_equal(unname(colMeans(xy)), c(0, 0))
})

test_that("coordinate problems give clear messages", {
  mk <- function(coords) SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(1, 1, nrow(coords))), spatialCoords = coords)
  expect_warning(tweedspot_model_coords(mk(cbind(x = 1:6, y = c(1, 3, 2, 5, 4, 6), z = 1))),
                 "uses the first two")
  expect_error(tweedspot_model_coords(mk(cbind(x = c(1:5, NA), y = 1:6))),
               "missing spatial coordinates")
  expect_error(tweedspot_model_coords(mk(cbind(x = 1:6, y = rep(2, 6)))), "is constant")
})

# ---- smaller input and output fixes -----------------------------------------

test_that("covariate names with spaces or brackets can be used", {
  counts <- storage_counts()
  spe <- storage_spe(counts)
  spe$`conc (mM)` <- rep(c(1, 2, 4), length.out = ncol(spe))
  spe$`batch id` <- rep(c("a", "b"), length.out = ncol(spe))
  cov <- tweedspot_covariates(spe, ~ log(`conc (mM)`) + `batch id`)
  expect_identical(cov$terms, "log(conc..mM.) + batch.id")
  expect_identical(names(cov$data), c("conc..mM.", "batch.id"))
  expect_identical(tweedspot_covariates(spe, c("conc (mM)", "batch id"))$terms,
                   "conc..mM. + batch.id")
  res <- tweedspot(spe, covariates = ~ `batch id` + `conc (mM)`, smooth_k = 8,
                   size_factor_method = "library_size",
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  expect_true(all(!is.na(SummarizedExperiment::rowData(res)$tweedspot_pval)))
  spe$x <- rep(1:2, length.out = ncol(spe))
  expect_error(tweedspot_covariates(spe, ~ x), "reserved names")
})

test_that("offset terms survive when a group drops a constant covariate", {
  cov <- list(
    data = data.frame(batch = factor(c("a", "a", "b", "b")), area = c(1, 2, 3, 4)),
    terms = "batch + offset(log(area))"
  )
  sub <- tweedspot_subset_covariates(cov, 1:2, "g")
  expect_identical(sub$terms, "offset(log(area))")
  expect_identical(names(sub$data), "area")
})

test_that("location masks with NA, empty results, and group selectors are handled", {
  counts <- storage_counts()
  spe <- storage_spe(counts)
  spe$keep <- c(NA, rep(c(1, 0), length.out = ncol(spe) - 1))
  res <- filter_features(spe, min_total = 0, min_pct_locations = 0, model_aware = FALSE,
                         location_mask = "keep")
  expect_identical(ncol(res), sum(spe$keep == 1, na.rm = TRUE))
  expect_warning(filter_features(spe, min_total = 1e9), "kept 0 feature")
})

test_that("plot files from cleaned names cannot collide", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  counts <- storage_counts()
  rownames(counts) <- c("A/B", "A_B", "flat")
  spe <- tweedspot(storage_spe(counts), smooth_k = 8, size_factor_method = "library_size",
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  dir <- tempfile("collide")
  f1 <- attr(plot_feature(spe, 1, output_dir = dir, formats = "png", dpi = 50), "files")
  f2 <- attr(plot_feature(spe, 2, output_dir = dir, formats = "png", dpi = 50), "files")
  expect_identical(basename(f1), "A_B_row1_tweedspot.png")
  expect_identical(basename(f2), "A_B_tweedspot.png")
  expect_message(plot_feature(spe, 2, output_dir = dir, formats = "png", dpi = 50),
                 "Overwriting existing file")
  unlink(dir, recursive = TRUE)
})

test_that("plot refits use the size factors and seed of the run", {
  counts <- storage_counts()
  spe <- tweedspot(storage_spe(counts), smooth_k = 8, size_factor_method = "library_size",
                   seed = 7, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  stored <- S4Vectors::metadata(spe)$tweedspot
  expect_identical(stored$seed, 7)
  # Size factors added after the run must not change the plot refit.
  spe$sizeFactor <- rep(c(0.5, 2), length.out = ncol(spe))
  settings <- tweedspot_plot_settings(spe)
  expect_identical(settings$size_factors, stored$size_factors)
  expect_identical(settings$seed, 7)
  # Overriding the normalization uses the new choice instead.
  expect_null(tweedspot_plot_settings(spe, list(size_factor_method = "scran"))$size_factors)
})

test_that("integer group labels can be selected as numbers", {
  set.seed(4)
  n <- 240
  counts <- rbind(a = stats::rpois(n, 5), b = stats::rpois(n, 3))
  colnames(counts) <- paste0("c", seq_len(n))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = cbind(x = stats::runif(n), y = stats::runif(n)),
    colData = S4Vectors::DataFrame(cluster = rep(1:2, each = n / 2))
  )
  res <- tweedspot(spe, group = "cluster", smooth_k = 8, min_group_size = 20,
                   size_factor_method = "library_size",
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  expect_identical(unique(get_tweedspot_results(res, group = 1)$group), "1")
  expect_identical(nrow(top_spatial_genes(res, n = 1, group = 2L)), 1L)
})
