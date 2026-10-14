spanorm_spe <- function(seed = 1) {
  set.seed(seed)
  grid <- expand.grid(x = seq_len(14), y = seq_len(14))
  n <- nrow(grid)
  tech <- exp(0.5 * as.numeric(scale(grid$x)))
  base <- exp(stats::rnorm(40, 1, 0.7))
  mu <- outer(base, tech)
  mu[1:5, grid$y > 9] <- 3 * mu[1:5, grid$y > 9]
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 5), nrow(mu))
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("spot", seq_len(n))
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(grid)
  )
}

# Fit SpaNorm once and cache it in the metadata; later calls reuse the fit.
spanorm_cached_spe <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      spe <- spanorm_spe()
      sn <- tweedspot_run_spanorm(spe, seed = 1)
      S4Vectors::metadata(spe)$SpaNorm <- sn$fit
      S4Vectors::metadata(spe)$SpaNorm_fingerprint <- sn$fingerprint
      cache <<- spe
    }
    cache
  }
})

test_that("PAC back-transform returns non-negative integers", {
  logpac <- log2(matrix(c(0, 1, 5, 12), 2) + 1)
  pac <- tweedspot_spanorm_pac(logpac)
  expect_equal(pac, matrix(c(0, 1, 5, 12), 2))
  expect_true(all(pac == round(pac)))
})

test_that("SpaNorm size factors are gene- and spot-specific and centred", {
  skip_if_not_installed("SpaNorm")
  spe <- spanorm_cached_spe()
  sn <- tweedspot_run_spanorm(spe, seed = 1)
  sf <- tweedspot_spanorm_offsets(sn$fit, rownames(spe), colnames(spe))
  expect_equal(dim(sf), dim(spe))
  expect_true(all(is.finite(sf) & sf > 0))
  expect_equal(unname(rowMeans(log(sf))), rep(0, nrow(spe)), tolerance = 1e-8)
  expect_null(sn$counts)
})

test_that("SpaNorm normalized counts are valid count responses", {
  skip_if_not_installed("SpaNorm")
  spe <- spanorm_cached_spe()
  fit_input <- tweedspot_prepare_input(
    spe, "counts", NULL, "spanorm", list(), spanorm_output = "counts", seed = 1
  )
  expect_equal(dim(fit_input$Y), dim(spe))
  expect_true(all(fit_input$Y >= 0 & fit_input$Y == round(fit_input$Y)))
  expect_equal(fit_input$libsz, rep(1, ncol(spe)))
  expect_false(is.null(fit_input$spanorm_fit))
})

test_that("tweedspot runs with both SpaNorm outputs and caches the fit", {
  skip_if_not_installed("SpaNorm")
  spe <- spanorm_cached_spe()
  for (out in c("size_factor", "counts")) {
    res <- tweedspot(
      spe,
      size_factor_method = "spanorm",
      spanorm_output = out,
      smooth_k = 10,
      seed = 1,
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    )
    rd <- SummarizedExperiment::rowData(res)
    expect_true(sum(!is.na(rd$tweedspot_pval)) >= 35)
    expect_false(is.null(S4Vectors::metadata(res)$SpaNorm))
  }
  # The cached fit is reused for the size-factor output without refitting.
  cached <- tweedspot_run_spanorm(res, seed = 1)
  expect_identical(cached$fit, S4Vectors::metadata(res)$SpaNorm)
})

test_that("a user-supplied adj.method is ignored with a warning", {
  skip_if_not_installed("SpaNorm")
  spe <- spanorm_cached_spe()
  expect_warning(
    tweedspot_run_spanorm(spe, spanorm_args = list(adj.method = "pearson"), seed = 1),
    "adj.method"
  )
})

test_that("a cached SpaNorm fit is reused only for the same data, settings and seed", {
  spe <- spanorm_spe()
  SummarizedExperiment::assay(spe, "other") <- SummarizedExperiment::assay(spe) + 1
  fp <- function(x, assay = "counts", args = list(), seed = 1) {
    tweedspot_spanorm_fingerprint(x, assay, args, seed)
  }
  base <- fp(spe)
  expect_identical(fp(spe), base)
  expect_false(identical(fp(spe[, rev(seq_len(ncol(spe)))]), base))
  expect_false(identical(fp(spe, assay = "other"), base))
  expect_false(identical(fp(spe, args = list(df.tps = 4)), base))
  expect_false(identical(fp(spe, seed = 2), base))
  # A fit without a TweedSpot fingerprint (or with another one) is not reused.
  skip_if_not_installed("SpaNorm")
  cached <- spanorm_cached_spe()
  stale <- cached
  S4Vectors::metadata(stale)$SpaNorm_fingerprint$seed <- 99
  fp_cached <- S4Vectors::metadata(cached)$SpaNorm_fingerprint
  expect_identical(fp_cached, fp(cached))
  expect_false(identical(S4Vectors::metadata(stale)$SpaNorm_fingerprint, fp(stale)))
})

test_that("the seed makes runs reproducible and restores the caller's RNG", {
  set.seed(5)
  before <- .Random.seed
  restore <- tweedspot_set_seed(1)
  a <- stats::runif(3)
  restore()
  expect_identical(.Random.seed, before)
  restore <- tweedspot_set_seed(1)
  b <- stats::runif(3)
  restore()
  expect_identical(a, b)
  expect_error(tweedspot_set_seed(1.5), "whole number")
  expect_error(tweedspot_set_seed("a"), "whole number")
})
