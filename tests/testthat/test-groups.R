# Group-specific (cell-type-specific) spatial testing.

group_spe <- function(n = 600, seed = 3) {
  set.seed(seed)
  x <- stats::runif(n, 0, 10)
  y <- stats::runif(n, 0, 10)
  # A on the left, B on the right, C scattered everywhere.
  type <- ifelse(stats::runif(n) < 0.25, "C", ifelse(x < 5, "A", "B"))
  sf <- exp(stats::rnorm(n, 0, 0.3))
  mk <- function(mu) stats::rnbinom(n, mu = mu * sf, size = 4)
  counts <- rbind(
    MARKER = mk(ifelse(type == "A", 12, 1)),               # composition only
    GRAD_A = mk(ifelse(type == "A", 2 + 8 * (y / 10), 4)), # gradient within A only
    SHARED = mk(2 + 8 * (y / 10)),                         # gradient in every type
    NULL1 = mk(3),
    NULL2 = mk(6)
  )
  colnames(counts) <- paste0("cell", seq_len(n))
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = cbind(x = x, y = y),
    colData = S4Vectors::DataFrame(
      cell_type = type,
      sizeFactor = sf,
      batch = ifelse(type == "C", "b1", sample(c("b1", "b2"), n, replace = TRUE))
    )
  )
}

run_groups <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- tweedspot(group_spe(), group = "cell_type", smooth_k = 10,
                          BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    }
    cache
  }
})

pair <- function(res, f, g) res[res$feature == f & res$group == g, ]

test_that("a pure composition effect is not called within any group", {
  overall <- tweedspot(group_spe(), smooth_k = 10,
                       BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  # The ordinary test sees the marker's spatial pattern...
  expect_lt(SummarizedExperiment::rowData(overall)["MARKER", "tweedspot_padj"], 1e-6)
  # ...but within every cell type it is flat.
  res <- get_tweedspot_results(run_groups(), group = TRUE)
  marker <- res[res$feature == "MARKER", ]
  expect_true(all(marker$tested))
  expect_true(all(marker$padj > 0.05))
})

test_that("within-group patterns are found only where they exist", {
  res <- get_tweedspot_results(run_groups(), group = TRUE)
  expect_lt(pair(res, "GRAD_A", "A")$padj, 1e-6)
  expect_gt(pair(res, "GRAD_A", "B")$padj, 0.05)
  expect_true(all(res$padj[res$feature == "SHARED"] < 1e-6))
  expect_true(all(res$padj[grepl("^NULL", res$feature)] > 0.05))
  rd <- SummarizedExperiment::rowData(run_groups())
  expect_match(rd["GRAD_A", "tweedspot_group_sig"], "(^|;)A(;|$)")
  expect_identical(rd["SHARED", "tweedspot_group_sig"], "A;B;C")
  expect_identical(rd["MARKER", "tweedspot_group_sig"], "")
})

test_that("group results have one row per feature and group and are stored", {
  spe <- run_groups()
  res <- S4Vectors::metadata(spe)$tweedspot_group$results
  expect_identical(nrow(res), 5L * 3L)
  expect_true(all(c("feature", "group", "n_locations", "n_detected", "tested", "reason",
                    "stat", "pval", "padj", "padj_feature", "edf", "dev_expl", "p",
                    "family") %in% names(res)))
  expect_identical(sort(unique(res$group)), c("A", "B", "C"))
  stored <- S4Vectors::metadata(spe)$tweedspot
  expect_identical(stored$group, "cell_type")
  expect_identical(stored$group_adjust, "stagewise")
  # Group mode leaves the ordinary per-feature columns untouched.
  expect_false("tweedspot_padj" %in% colnames(SummarizedExperiment::rowData(spe)))
})

test_that("small groups and unsupported pairs are recorded, not tested", {
  spe <- group_spe()
  counts <- SummarizedExperiment::assay(spe)
  type <- spe$cell_type
  # NULL1 is detected in only 5 cells of type B.
  b <- which(type == "B")
  counts["NULL1", b] <- 0
  counts["NULL1", b[1:5]] <- 3
  SummarizedExperiment::assay(spe) <- counts
  res <- get_tweedspot_results(
    tweedspot(spe, group = "cell_type", smooth_k = 10, min_group_size = 160,
              BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
    group = TRUE
  )
  small <- names(which(table(type) < 160))
  expect_true(length(small) >= 1L)
  expect_true(all(!res$tested[res$group %in% small]))
  expect_true(all(res$reason[res$group %in% small] == "group_too_small"))
  b_row <- pair(res, "NULL1", "B")
  if (!"B" %in% small) {
    expect_false(b_row$tested)
    expect_identical(b_row$reason, "low_model_support")
  }
  expect_true(all(is.na(res$padj[!res$tested])))
})

test_that("covariates constant within a group are dropped for that group", {
  expect_message(
    tweedspot(group_spe(), group = "cell_type", covariates = ~ batch, smooth_k = 10,
              BPPARAM = BiocParallel::SerialParam(), verbose = TRUE),
    "Group 'C': dropping covariate"
  )
  cov <- list(
    data = data.frame(batch = factor(c("a", "a", "b", "b")), age = c(1, 2, 3, 4)),
    terms = "batch + age + batch:age"
  )
  sub <- tweedspot_subset_covariates(cov, 1:2, "g")
  expect_identical(sub$terms, "age")
  expect_identical(names(sub$data), "age")
  expect_identical(tweedspot_subset_covariates(cov, 1:4, "g")$terms, cov$terms)
  expect_null(tweedspot_subset_covariates(list(data = cov$data[, "batch", drop = FALSE],
                                               terms = "batch"), 1:2, "g"))
  expect_null(tweedspot_subset_covariates(NULL, 1:2, "g"))
})

test_that("group p-value adjustments behave as documented", {
  res <- data.frame(
    feature = rep(c("f1", "f2", "f3"), each = 2),
    feature_index = rep(1:3, each = 2),
    group = rep(c("A", "B"), 3),
    tested = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    pval = c(1e-6, 0.4, 0.03, 0.04, 0.5, NA)
  )
  ok <- res$tested
  global <- tweedspot_group_adjust(res, "global", "BH")
  expect_equal(global$padj[ok], stats::p.adjust(res$pval[ok], "BH"))
  per <- tweedspot_group_adjust(res, "per_group", "BH")
  expect_equal(per$padj[per$group == "A"], stats::p.adjust(res$pval[res$group == "A"], "BH"))

  stage <- tweedspot_group_adjust(res, "stagewise", "BH", alpha = 0.05)
  screen <- stats::p.adjust(c(cct(c(1e-6, 0.4)), cct(c(0.03, 0.04)), 0.5), "BH")
  expect_equal(stage$padj_feature[stage$feature == "f1"], rep(screen[1], 2))
  # Confirmation is at alpha * R / m: here R = 1 of m = 3 features pass.
  expect_identical(sum(screen <= 0.05), 1L)
  holm_f1 <- pmin(1, stats::p.adjust(c(1e-6, 0.4), "holm") * 3 / 1)
  expect_equal(stage$padj[1:2], pmax(holm_f1, screen[1]))
  expect_true(all(stage$padj[ok] >= stage$pval[ok]))
  expect_true(is.na(stage$padj[6]))
})

test_that("stage-wise adjustment controls the pair-level false discovery rate", {
  set.seed(10)
  nf <- 1500
  G <- 4
  true_f <- sample(nf, 150)
  true_g <- sample(G, 150, replace = TRUE)
  fdp <- vapply(1:15, function(rep) {
    p <- matrix(stats::runif(nf * G), nf, G)
    p[cbind(true_f, true_g)] <- stats::rbeta(150, 0.05, 1)
    res <- data.frame(feature = rep(paste0("f", seq_len(nf)), G),
                      feature_index = rep(seq_len(nf), G),
                      group = rep(LETTERS[seq_len(G)], each = nf),
                      tested = TRUE, pval = as.vector(p))
    truth <- rep(FALSE, nf * G)
    truth[(true_g - 1) * nf + true_f] <- TRUE
    adj <- tweedspot_group_adjust(res, "stagewise", "BH", alpha = 0.05)
    hit <- adj$padj <= 0.05
    if (any(hit)) sum(hit & !truth) / sum(hit) else 0
  }, numeric(1))
  expect_lt(mean(fdp), 0.065)
})

test_that("duplicated feature names stay separate in group results", {
  spe <- group_spe()
  rownames(spe)[5] <- rownames(spe)[4]
  res <- tweedspot(spe, group = "cell_type", smooth_k = 10,
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  tab <- get_tweedspot_results(res, group = TRUE, sort_by = "none")
  expect_identical(sort(unique(tab$feature_index)), 1:5)
  rd <- SummarizedExperiment::rowData(res)
  expect_length(rd$tweedspot_group_padj, 5L)
})

test_that("group mode works without row names", {
  spe <- group_spe()
  rownames(spe) <- NULL
  res <- tweedspot(spe, group = "cell_type", smooth_k = 10,
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  tab <- get_tweedspot_results(res, group = TRUE)
  expect_identical(nrow(tab), 15L)
  expect_true(all(grepl("^feature[0-9]+$", tab$feature)))
})

test_that("accessors return group results", {
  spe <- run_groups()
  a <- get_tweedspot_results(spe, group = "A")
  expect_identical(unique(a$group), "A")
  expect_identical(nrow(a), 5L)
  expect_true(is.unsorted(a$padj, na.rm = TRUE) == FALSE)
  ab <- get_tweedspot_results(spe, group = c("A", "B"), sort_by = "none")
  expect_identical(sort(unique(ab$group)), c("A", "B"))
  expect_error(get_tweedspot_results(spe, group = "Z"), "Unknown group")
  expect_error(get_tweedspot_results(spe, group = 1), "group")
  top <- top_spatial_genes(spe, n = 1, group = c("A", "C"))
  expect_identical(top$group, c("A", "C"))
  expect_identical(top$feature[top$group == "A"], top_spatial_genes(spe, n = 1, group = "A")$feature)
  expect_identical(nrow(top_spatial_genes(spe, n = 4, group = TRUE)), 4L)
  plain <- tweedspot(group_spe(), smooth_k = 10, BPPARAM = BiocParallel::SerialParam(),
                     verbose = FALSE)
  expect_error(get_tweedspot_results(plain, group = TRUE), "run `tweedspot\\(group")
})

test_that("group labels are validated and missing labels are skipped", {
  spe <- group_spe()
  expect_error(tweedspot(spe, group = "missing", verbose = FALSE), "must name a column")
  spe$score <- stats::runif(ncol(spe))
  expect_error(tweedspot(spe, group = "score", verbose = FALSE), "discrete labels")
  spe$cell_type[1:20] <- NA
  expect_message(
    res <- tweedspot(spe, group = "cell_type", smooth_k = 10,
                     BPPARAM = BiocParallel::SerialParam(), verbose = TRUE),
    "20 location\\(s\\) with a missing"
  )
  tab <- get_tweedspot_results(res, group = TRUE)
  expect_identical(sum(tab$n_locations[tab$feature == "NULL1"]), ncol(spe) - 20L)
})

test_that("a failing group gives a warning, not an error", {
  calls <- 0L
  real <- tweedspot_agnostic
  local_mocked_bindings(tweedspot_agnostic = function(Y, coords, ...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("simulated failure")
    real(Y, coords, ...)
  })
  expect_warning(
    res <- tweedspot(group_spe(), group = "cell_type", smooth_k = 10,
                     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
    "All TweedSpot fits failed in group 'B'"
  )
  tab <- get_tweedspot_results(res, group = TRUE)
  expect_true(all(tab$reason[tab$group == "B"] == "fit_failed"))
  expect_true(all(tab$tested[tab$group == "A"]))
})

test_that("filter_features checks model support within each group", {
  spe <- group_spe()
  counts <- SummarizedExperiment::assay(spe)
  type <- spe$cell_type
  a <- which(type == "A")
  # ONLY_A: detected in 30 cells, all of type A.
  only_a <- rep(0, ncol(spe))
  only_a[a[1:30]] <- 5
  # NOWHERE: detected in 4 cells of each type.
  nowhere <- rep(0, ncol(spe))
  nowhere[unlist(lapply(c("A", "B", "C"), function(g) which(type == g)[1:4]))] <- 5
  counts <- rbind(counts, ONLY_A = only_a, NOWHERE = nowhere)
  spe2 <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = SpatialExperiment::spatialCoords(spe),
    colData = SummarizedExperiment::colData(spe)
  )
  base <- list(input = spe2, min_total = 0, min_pct_locations = 0, smooth_k = 10)
  grouped <- do.call(filter_features, c(base, group = "cell_type", action = "flag"))
  rd <- SummarizedExperiment::rowData(grouped)
  expect_true(rd["ONLY_A", "filter_keep"])
  expect_false(rd["NOWHERE", "filter_keep"])
  expect_identical(rd["NOWHERE", "filter_reason"], "low_model_support")
  support <- S4Vectors::metadata(grouped)$filter_features$group_support
  expect_identical(colnames(support), c("A", "B", "C"))
  expect_identical(unname(support["ONLY_A", ]), c(TRUE, FALSE, FALSE))
  # Globally NOWHERE's 12 detections pass a 10-function smooth; within each
  # group (4 detections each) they do not.
  expect_true("NOWHERE" %in% rownames(do.call(filter_features, base)))
  # Groups below min_group_size are ignored.
  big_only <- do.call(filter_features, c(base, group = "cell_type", min_group_size = 160,
                                         action = "flag"))
  expect_false("C" %in% colnames(S4Vectors::metadata(big_only)$filter_features$group_support))
})

# ---- results from different runs never mix -----------------------------------

test_that("a group run replaces an earlier ordinary run, and vice versa", {
  spe <- group_spe()
  ordinary <- tweedspot(spe, smooth_k = 10, BPPARAM = BiocParallel::SerialParam(),
                        verbose = FALSE)
  expect_message(
    grouped <- tweedspot(ordinary, group = "cell_type", smooth_k = 10,
                         BPPARAM = BiocParallel::SerialParam(), verbose = TRUE),
    "Replacing results from a previous tweedspot\\(\\) run"
  )
  rd <- SummarizedExperiment::rowData(grouped)
  expect_false(any(c("tweedspot_padj", "tweedspot_pval", "tweedspot_stat") %in% colnames(rd)))
  expect_error(get_tweedspot_results(grouped), "result columns not found")

  back <- tweedspot(grouped, smooth_k = 10, BPPARAM = BiocParallel::SerialParam(),
                    verbose = FALSE)
  expect_null(S4Vectors::metadata(back)$tweedspot_group)
  expect_false(any(grepl("^tweedspot_group", colnames(SummarizedExperiment::rowData(back)))))
  expect_error(get_tweedspot_results(back, group = TRUE), "run `tweedspot\\(group")
  expect_null(S4Vectors::metadata(back)$tweedspot$group)
})

test_that("group results follow feature subsetting and warn on location changes", {
  spe <- run_groups()
  sub <- spe[c("SHARED", "MARKER"), ]
  tab <- get_tweedspot_results(sub, group = TRUE, sort_by = "none")
  expect_identical(sort(unique(tab$feature)), c("MARKER", "SHARED"))
  expect_identical(nrow(tab), 6L)
  expect_identical(tab$feature_index, match(tab$feature, rownames(sub)))
  expect_warning(get_tweedspot_results(spe[, 1:100], group = TRUE),
                 "locations were subset or reordered")
  ordinary <- tweedspot(group_spe(), smooth_k = 10,
                        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  expect_warning(get_tweedspot_results(ordinary[, 1:100]), "locations were subset")
  expect_silent(get_tweedspot_results(ordinary[1:2, ]))
  expect_identical(get_tweedspot_results(ordinary, sort_by = "none")$feature,
                   rownames(ordinary))
})

test_that("plots follow a group run and refuse mismatched requests", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spe <- run_groups()
  expect_error(plot_feature(spe, "GRAD_A"), "Choose the group to plot")
  expect_error(plot_feature(spe, "GRAD_A", group = "Z"), "must be one of")
  p <- plot_feature(spe, "GRAD_A", group = "A")
  expect_identical(p$patches$annotation$title, "GRAD_A in A")
  # Only cells of type A are drawn.
  expect_identical(nrow(p[[1]]$data), sum(spe$cell_type == "A"))

  top <- plot_top_features(spe, n_top = 2, group = "A", label_stats = TRUE)
  expect_identical(top$labels$title, "A")
  expected <- get_tweedspot_results(spe, group = "A")
  expected <- expected[expected$tested, ]
  expect_match(levels(top$data$panel)[1], paste0("^", expected$feature[1], "\\n"))
  expect_identical(nrow(top$data), 2L * sum(spe$cell_type == "A"))

  ordinary <- tweedspot(group_spe(), smooth_k = 10,
                        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  expect_error(plot_feature(ordinary, "GRAD_A", group = "A"), "applies only to results")
  expect_error(plot_top_features(spe), "Choose the group to plot")
})

test_that("automatic point size leaves the random-number stream untouched", {
  set.seed(99)
  before <- .Random.seed
  tweedspot_auto_point_size(as.matrix(expand.grid(x = 1:50, y = 1:50)), 2, 2)
  expect_identical(.Random.seed, before)
})
