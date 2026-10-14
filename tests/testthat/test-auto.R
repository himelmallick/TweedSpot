# Tests for `tweedie_p = "auto"`: per-gene choice of the Tweedie variance power
# by an information criterion.

auto_frame <- function(family, seed = 1) {
  set.seed(seed)
  grid <- expand.grid(x = seq(-1, 1, length.out = 20), y = seq(-1, 1, length.out = 20))
  n <- nrow(grid)
  mu <- exp(log(10) + 0.6 * sin(2 * grid$x) + 0.4 * grid$y)
  expr <- switch(family,
    poisson = stats::rpois(n, mu),
    nb = stats::rnbinom(n, mu = mu, size = 1.5),
    gaussian = stats::rnorm(n, mu + 40, 2),
    gamma = stats::rgamma(n, shape = 2, rate = 2 / mu),
    inverse_gaussian = {
      # Inverse Gaussian draws (Michael, Schucany & Haas, 1976), dispersion 0.05.
      lambda <- 1 / 0.05
      nu <- stats::rnorm(n)^2
      x <- mu + mu^2 * nu / (2 * lambda) -
        mu / (2 * lambda) * sqrt(4 * mu * lambda * nu + mu^2 * nu^2)
      ifelse(stats::runif(n) <= mu / (mu + x), x, mu^2 / x)
    },
    tweedie = {
      # Compound Poisson-gamma with p = 1.5, phi = 2.
      p <- 1.5
      phi <- 2
      lam <- mu^(2 - p) / (phi * (2 - p))
      shape <- (2 - p) / (p - 1)
      scale <- phi * (p - 1) * mu^(p - 1)
      N <- stats::rpois(n, lam)
      vapply(seq_len(n), function(i) {
        if (N[i] == 0) 0 else sum(stats::rgamma(N[i], shape = shape, scale = scale[i]))
      }, numeric(1))
    }
  )
  data.frame(expr = expr, x = grid$x, y = grid$y, libsz = 1)
}

auto_spe <- function(counts) {
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
  side <- sqrt(ncol(counts))
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  )
}

auto_formula <- expr ~ s(x, y, k = 20) + offset(log(libsz))

fit_auto <- function(d, criterion = "BIC") {
  tweedspot_fit_auto(auto_formula, d, "REML", FALSE, TRUE, 1L, criterion)
}

test_that("tweedie_p = 'auto' is recognised and validated", {
  expect_identical(tweedspot_resolve_family("tw", "auto"), "auto")
  expect_identical(tweedspot_resolve_family("tw", "AUTO"), "auto")
  expect_error(tweedspot_resolve_family("tw", "bic"), "auto")
  expect_error(tweedspot_resolve_family("tw", c("auto", "auto")), "auto")
})

test_that("fixed variance powers map to log-link families", {
  y <- c(0, 1, 2, 3)
  expect_identical(tweedspot_resolve_family("tw", 0)$family, "gaussian")
  expect_identical(tweedspot_resolve_family("tw", 0)$link, "log")
  expect_identical(tweedspot_resolve_family("tw", 1)$family, "poisson")
  expect_match(tweedspot_resolve_family("tw", 1.5)$family, "^Tweedie")
  expect_error(tweedspot_resolve_family("tw", 2, Y = y), "strictly positive")
  expect_error(tweedspot_resolve_family("tw", 3, Y = y), "strictly positive")
  expect_identical(tweedspot_resolve_family("tw", 3, Y = y + 1)$family, "inverse.gaussian")
})

test_that("candidate set depends on zeros and integer support", {
  expect_named(tweedspot_auto_candidates(c(0, 1, 4, 2)),
               c("gaussian", "tweedie"))
  expect_named(tweedspot_auto_candidates(c(3, 1, 4, 2)),
               c("gaussian", "tweedie", "gamma", "inverse_gaussian"))
  expect_named(tweedspot_auto_candidates(c(0, 1.5, 4.2)),
               c("gaussian", "tweedie"))
  expect_named(tweedspot_auto_candidates(c(0.3, 1.5, 4.2)),
               c("gaussian", "tweedie", "gamma", "inverse_gaussian"))
  expect_true(tweedspot_is_integer_valued(c(0, 2, NA, 5)))
  expect_false(tweedspot_is_integer_valued(c(0, 2.5)))
})

test_that("every candidate shares the log link", {
  cands <- tweedspot_auto_candidates(c(3, 1, 4, 2))
  expect_true(all(vapply(cands, function(f) f$link, character(1)) == "log"))
})

test_that("information criteria match mgcv's log-likelihood and df", {
  d <- auto_frame("nb")
  fit <- mgcv::gam(auto_formula, family = tweedspot_tw(), data = d, method = "REML")
  ll <- stats::logLik(fit)
  expect_equal(
    tweedspot_information_criterion(fit, "BIC"),
    -2 * as.numeric(ll) + log(nrow(d)) * attr(ll, "df")
  )
  expect_equal(tweedspot_information_criterion(fit, "AIC"), stats::AIC(fit))
  # The estimated Tweedie power is counted as an extra parameter.
  fit_fixed <- mgcv::gam(auto_formula, family = mgcv::Tweedie(1.5, link = "log"),
                         data = d, method = "REML")
  expect_gt(attr(ll, "df") - sum(fit$edf), attr(stats::logLik(fit_fixed), "df") - sum(fit_fixed$edf))
})

test_that("auto keeps the candidate with the smallest criterion", {
  d <- auto_frame("nb")
  fit <- fit_auto(d)
  ic <- attr(fit, "tweedspot_ic")
  expect_named(ic, names(tweedspot_auto_candidates(d$expr)))
  expect_identical(attr(fit, "tweedspot_selected"), names(ic)[which.min(ic)])
})

test_that("auto recovers the generating family", {
  expected <- c(
    poisson = "tweedie",
    nb = "tweedie",
    gaussian = "gaussian",
    gamma = "gamma",
    inverse_gaussian = "inverse_gaussian",
    tweedie = "tweedie"
  )
  for (fam in names(expected)) {
    for (seed in 1:2) {
      sel <- attr(fit_auto(auto_frame(fam, seed)), "tweedspot_selected")
      expect_identical(sel, expected[[fam]], info = paste(fam, "seed", seed))
    }
  }
})

test_that("auto never selects Poisson and keeps Tweedie at p >= 1.1 on counts", {
  # Poisson counts are represented by Tweedie at the lower bound. Without the
  # p >= 1.1 bound, Tweedie with p = 1.01 spikes at the integers (seeds 2, 4,
  # 5 and 6), which inflates its likelihood and the spatial test.
  for (seed in c(2, 4, 5, 6)) {
    fit <- fit_auto(auto_frame("poisson", seed))
    expect_false("poisson" %in% names(attr(fit, "tweedspot_ic")))
    expect_identical(attr(fit, "tweedspot_selected"), "tweedie", info = paste("seed", seed))
    expect_gte(tweedspot_family_power(fit), tweedspot_auto_min_p() - 1e-8)
  }
})

test_that("auto is calibrated on null counts with little overdispersion", {
  # These are the genes on which a Poisson candidate used to be selected and
  # then falsely called spatially variable.
  set.seed(31)
  grid <- expand.grid(x = seq(-1, 1, length.out = 20), y = seq(-1, 1, length.out = 20))
  pvals <- vapply(1:60, function(i) {
    d <- data.frame(expr = stats::rnbinom(400, mu = stats::runif(1, 4, 8),
                                          size = stats::runif(1, 20, 40)),
                    x = grid$x, y = grid$y, libsz = 1)
    summary(fit_auto(d))$s.table[1, "p-value"]
  }, numeric(1))
  expect_lte(mean(pvals < 0.05), 0.12)
})

test_that("estimated Tweedie power is reported and close to the truth", {
  fit <- fit_auto(auto_frame("tweedie"))
  expect_identical(attr(fit, "tweedspot_selected"), "tweedie")
  expect_equal(tweedspot_family_power(fit), 1.5, tolerance = 0.1)
  expect_gte(tweedspot_family_power(fit_auto(auto_frame("nb"))), tweedspot_auto_min_p())
})

test_that("tweedspot_family_power maps each family to its power", {
  d <- auto_frame("gamma")
  pw <- function(fam) {
    tweedspot_family_power(mgcv::gam(auto_formula, family = fam, data = d, method = "REML"))
  }
  expect_identical(pw(stats::gaussian(link = "log")), 0)
  expect_identical(pw(stats::Gamma(link = "log")), 2)
  expect_identical(pw(stats::inverse.gaussian(link = "log")), 3)
  expect_equal(pw(mgcv::Tweedie(1.3, link = "log")), 1.3)
  d$expr <- round(d$expr)
  expect_identical(pw(stats::poisson(link = "log")), 1)
  expect_identical(tweedspot_power_label(c(0, 1, 1.4, 2, 3, NA)),
                   c("gaussian", "poisson", "tweedie", "gamma", "inverse_gaussian", NA))
})

test_that("tweedspot() with tweedie_p = 'auto' reports the chosen power per gene", {
  set.seed(2)
  spe <- auto_spe(matrix(c(
    stats::rpois(36, 5),
    stats::rnbinom(36, mu = 5, size = 1),
    stats::rpois(36, 20)
  ), nrow = 3, byrow = TRUE))
  for (crit in c("BIC", "AIC")) {
    res <- tweedspot(
      spe,
      tweedie_p = "auto",
      p_criterion = crit,
      smooth_k = 5,
      size_factor_method = "library_size",
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    )
    rd <- SummarizedExperiment::rowData(res)
    expect_true(all(c("tweedspot_p", "tweedspot_family") %in% colnames(rd)))
    ok <- !is.na(rd$tweedspot_p)
    expect_true(any(ok))
    expect_true(all(rd$tweedspot_family[ok] %in%
      c("gaussian", "tweedie", "gamma", "inverse_gaussian")))
    expect_identical(rd$tweedspot_family[ok], tweedspot_power_label(rd$tweedspot_p[ok]))
  }
  expect_error(
    tweedspot(spe, tweedie_p = "auto", p_criterion = "DIC", verbose = FALSE),
    "should be one of"
  )
})

test_that("auto works with the two-part model and with plot_feature()", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(3)
  spe <- auto_spe(matrix(stats::rnbinom(4 * 36, mu = 3, size = 2), nrow = 4))
  res <- tweedspot(
    spe,
    tweedie_p = "auto",
    two_part = TRUE,
    smooth_k = 5,
    size_factor_method = "library_size",
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )
  expect_true(any(!is.na(SummarizedExperiment::rowData(res)$tweedspot_pval)))
  p <- plot_feature(spe, 1, tweedie_p = "auto", smooth_k = 5,
                    size_factor_method = "library_size")
  expect_s3_class(p, "patchwork")
})

test_that("auto passes on warnings from the selected fit only", {
  real_fit <- tweedspot_fit
  local_mocked_bindings(tweedspot_fit = function(formula, family, ...) {
    if (inherits(family, "family") && identical(family$family, "Gamma")) {
      warning("gamma candidate warned")
    }
    real_fit(formula, family, ...)
  })
  # Gamma data: Gamma is selected, so its warning is shown.
  expect_warning(fit_auto(auto_frame("gamma")), "gamma candidate warned")
  # Positive Gaussian data: Gamma is a candidate but Gaussian is selected, so
  # the Gamma candidate's warning is not shown.
  d <- auto_frame("gaussian")
  expect_true("gamma" %in% names(tweedspot_auto_candidates(d$expr)))
  expect_no_warning(sel <- fit_auto(d))
  expect_identical(attr(sel, "tweedspot_selected"), "gaussian")
})
