# Regression tests: a family object reused across genes must not carry one
# gene's estimated Tweedie power (or NB theta) into the next gene's fit.

state_data <- function() {
  set.seed(21)
  grid <- expand.grid(x = seq(-1, 1, length.out = 25), y = seq(-1, 1, length.out = 25))
  n <- nrow(grid)
  list(
    grid = grid,
    # A gene that drives the power to its lower bound, then an ordinary NB gene.
    first = stats::rpois(n, 6),
    second = stats::rnbinom(n, mu = 4, size = 2.5)
  )
}

fit_power <- function(expr, family, grid) {
  d <- data.frame(expr = expr, x = grid$x, y = grid$y, libsz = 1)
  fit <- tweedspot_fit(expr ~ s(x, y, k = 20) + offset(log(libsz)), family, d,
                       "REML", FALSE, TRUE, 1L)
  c(p = tweedspot_family_power(fit), pval = summary(fit)$s.table[1, "p-value"])
}

test_that("a reused Tweedie family gives order-independent fits", {
  sd <- state_data()
  shared <- tweedspot_resolve_family("tw", NULL, Y = sd$second)
  fit_power(sd$first, shared, sd$grid)
  after_first <- fit_power(sd$second, shared, sd$grid)
  alone <- fit_power(sd$second, tweedspot_resolve_family("tw", NULL, Y = sd$second), sd$grid)
  expect_equal(after_first, alone, tolerance = 1e-6)
})

test_that("a reused negative binomial family gives order-independent fits", {
  sd <- state_data()
  shared <- tweedspot_resolve_family(mgcv::nb(), NULL)
  fit_power(sd$first, shared, sd$grid)
  d <- data.frame(expr = sd$second, x = sd$grid$x, y = sd$grid$y, libsz = 1)
  form <- expr ~ s(x, y, k = 20) + offset(log(libsz))
  after_first <- tweedspot_fit(form, shared, d, "REML", FALSE, TRUE, 1L)
  alone <- tweedspot_fit(form, tweedspot_resolve_family(mgcv::nb(), NULL), d,
                         "REML", FALSE, TRUE, 1L)
  expect_equal(after_first$family$getTheta(TRUE), alone$family$getTheta(TRUE),
               tolerance = 1e-6)
})

test_that("tweedspot() results do not depend on gene order", {
  sd <- state_data()
  counts <- rbind(first = sd$first, second = sd$second)
  make_spe <- function(m) {
    colnames(m) <- paste0("loc", seq_len(ncol(m)))
    SpatialExperiment::SpatialExperiment(
      assays = list(counts = m),
      spatialCoords = as.matrix(sd$grid),
      colData = S4Vectors::DataFrame(sizeFactor = rep(1, ncol(m)))
    )
  }
  run <- function(m) {
    rd <- SummarizedExperiment::rowData(tweedspot(
      make_spe(m), smooth_k = 20, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
    ))
    stats::setNames(rd$tweedspot_pval, rownames(m))
  }
  forward <- run(counts)
  reverse <- run(counts[2:1, ])
  expect_equal(forward[c("first", "second")], reverse[c("first", "second")], tolerance = 1e-6)
})

test_that("default Tweedie power is bounded at 1.1 for integer data only", {
  bound <- function(fam) get(".a", envir = environment(fam$getTheta))
  expect_equal(bound(tweedspot_resolve_family("tw", NULL, Y = c(0, 1, 5))), 1.1)
  expect_equal(bound(tweedspot_resolve_family("tw", NULL, Y = c(0, 1.5, 5))), 1.01)
})
