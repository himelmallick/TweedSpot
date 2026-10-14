test_that("filter_features filters features and empty locations", {
  spe <- mock_spe()

  filtered <- filter_features(
    input = spe,
    assay_name = "counts",
    min_total = 8,
    min_pct_locations = 40,
    min_locations = 3,
    min_mean = 1,
    drop_empty_locations = TRUE
  )

  expect_s4_class(filtered, "SpatialExperiment")
  expect_true(nrow(filtered) < nrow(spe))
  expect_true(ncol(filtered) < ncol(spe))
  expect_true(all(colSums(SummarizedExperiment::assay(filtered, "counts")) > 0))
})

test_that("each threshold removes exactly the features that fail it", {
  counts <- rbind(
    low_total = c(1, 0, 0, 0, 1, 0),
    rare      = c(0, 0, 0, 0, 0, 9),
    flat      = c(2, 2, 2, 2, 2, 2),
    good      = c(3, 0, 5, 1, 4, 2)
  )
  spe <- mock_spe(counts)
  f <- function(...) {
    args <- utils::modifyList(
      list(input = spe, min_total = 0, min_pct_locations = 0, model_aware = FALSE),
      list(...)
    )
    rownames(do.call(filter_features, args))
  }

  expect_false("gene1" %in% f(min_total = 5))
  expect_false("gene2" %in% f(min_locations = 2))
  expect_false("gene2" %in% f(min_pct_locations = 30))
  expect_false("gene3" %in% f(min_var = 0.1))
  expect_false("gene1" %in% f(min_mean = 0.5))
  expect_identical(f(min_total = 5, min_locations = 2, min_var = 0.1), "gene4")
})

test_that("detection_threshold supports intensity data with a background floor", {
  intens <- rbind(
    background = c(0.2, 0.3, 0.1, 0.2, 0.3, 0.2),
    signal     = c(0.2, 5.1, 0.1, 6.0, 4.2, 0.3)
  )
  spe <- mock_spe(intens)
  kept <- filter_features(spe, min_total = 0, min_locations = 3,
                          detection_threshold = 1, model_aware = FALSE)
  expect_identical(rownames(kept), "gene2")
  expect_identical(
    nrow(filter_features(spe, min_total = 0, min_locations = 3, model_aware = FALSE)),
    2L
  )
})

test_that("filter_features rejects negative assays and bad arguments", {
  spe <- mock_spe(matrix(c(-1, 2, 3, 4, 5, 6), nrow = 1))
  expect_error(filter_features(spe), "non-negative")
  spe <- mock_spe()
  expect_error(filter_features(spe, min_pct_locations = 120), "min_pct_locations")
  expect_error(filter_features(spe, detection_threshold = -1), "detection_threshold")
  expect_error(filter_features(spe, exclude = 1:2), "exclude")
})

test_that("explicit exclusions are modality agnostic", {
  counts <- matrix(5, nrow = 4, ncol = 6)
  spe <- mock_spe(counts)
  rownames(spe) <- c("CD3E", "NegPrb1", "BLANK_0004", "Panck")
  SummarizedExperiment::rowData(spe)$label <- c("a", "b", "c", "DAPI")

  expect_identical(
    rownames(filter_features(spe, exclude = "CD3E")),
    c("NegPrb1", "BLANK_0004", "Panck")
  )
  expect_identical(
    rownames(filter_features(spe, exclude_pattern = c("^NegPrb", "^BLANK"))),
    c("CD3E", "Panck")
  )
  # Matching also uses feature_name_col when supplied.
  expect_identical(
    rownames(filter_features(spe, exclude = "DAPI", feature_name_col = "label")),
    c("CD3E", "NegPrb1", "BLANK_0004")
  )
})

test_that("filter_features reports retained and removed features", {
  spe <- mock_spe()

  expect_message(
    filter_features(
      input = spe,
      assay_name = "counts",
      min_total = 3,
      exclude_pattern = "^gene1$",
      report = TRUE
    ),
    "filter_features retained"
  )
  expect_message(
    filter_features(spe, min_total = 3, exclude = "gene1", exclude_mito = TRUE,
                    feature_name_col = "gene_name", report = TRUE),
    "excluded=1, mito_ribo=1"
  )
  expect_message(
    filter_features(spe, action = "flag", report = TRUE),
    "flagged to keep"
  )
})

test_that("filter_features can exclude mitochondrial and ribosomal genes", {
  spe <- mock_spe()

  filtered <- filter_features(
    input = spe,
    assay_name = "counts",
    exclude_mito = TRUE,
    exclude_ribo = TRUE,
    feature_name_col = "gene_name"
  )

  expect_false(any(rownames(filtered) %in% c("gene1", "gene2")))
})

test_that("filter_genes is no longer exported", {
  expect_false("filter_genes" %in% getNamespaceExports("TweedSpot"))
})

# ---- helpers for model-aware, flag, and location QC tests -------------------

grid_spe <- function(counts, side = 6) {
  rownames(counts) <- paste0("f", seq_len(nrow(counts)))
  colnames(counts) <- paste0("loc", seq_len(ncol(counts)))
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  )
}

sparse_counts <- function(n_detected, n = 36, value = 10) {
  x <- rep(0, n)
  x[seq_len(n_detected)] <- value
  x
}

test_that("model-aware minimum follows the smooth's basis dimension", {
  # Features detected in 3, 8, 12 and all 36 of 36 locations.
  spe <- grid_spe(rbind(sparse_counts(3), sparse_counts(8), sparse_counts(12),
                        rep(4, 36)))
  base <- list(input = spe, min_total = 0, min_pct_locations = 0)
  kept <- function(...) rownames(do.call(filter_features, c(base, list(...))))

  expect_identical(kept(smooth_k = 10), c("f3", "f4"))
  expect_identical(kept(smooth_k = 5), c("f2", "f3", "f4"))
  expect_identical(kept(smooth_k = 10, model_aware = FALSE), c("f1", "f2", "f3", "f4"))
  # Default basis is min(30, unique locations - 1) = 30 here.
  expect_identical(kept(), "f4")
  res <- filter_features(spe, min_total = 0, min_pct_locations = 0, smooth_k = 10)
  expect_identical(S4Vectors::metadata(res)$filter_features$min_support, 10L)
  expect_error(filter_features(spe, smooth_k = 2), "smooth_k")
})

test_that("model-aware support is counted on locations that pass QC", {
  counts <- rbind(sparse_counts(12), rep(4, 36))
  spe <- grid_spe(counts)
  SummarizedExperiment::colData(spe)$in_tissue <- c(rep(FALSE, 6), rep(TRUE, 30))
  res <- filter_features(spe, min_total = 0, min_pct_locations = 0, smooth_k = 10,
                         location_mask = "in_tissue")
  # f1 is detected in 12 locations, but only 6 of them are in tissue.
  expect_identical(rownames(res), "f2")
  expect_identical(ncol(res), 30L)
})

test_that("two-part reporting counts features whose binary part is skipped", {
  spe <- grid_spe(rbind(sparse_counts(20), rep(4, 36)))
  expect_message(
    res <- filter_features(spe, min_total = 0, min_pct_locations = 0, smooth_k = 10,
                           two_part = TRUE, report = TRUE),
    "1 kept feature"
  )
  expect_identical(nrow(res), 2L)
  expect_identical(S4Vectors::metadata(res)$filter_features$binary_part_skipped, 1L)
})

test_that("binary-part support rule and NA-safe combination", {
  expect_true(tweedspot_binary_estimable(sparse_counts(20), 10))
  expect_false(tweedspot_binary_estimable(rep(4, 36), 10))
  expect_false(tweedspot_binary_estimable(sparse_counts(5), 10))
  M <- cbind(pval = c(0.01, 0.2, NA), plogit = c(NA, 0.3, 0.4))
  combined <- tweedspot_combine_pvalues(M)
  expect_equal(combined[1], 0.01, tolerance = 1e-8)
  expect_equal(combined[2], cct(c(0.2, 0.3)))
  expect_true(is.na(combined[3]))
})

test_that("two-part tweedspot uses the count part alone when zeros are too few", {
  set.seed(4)
  grid <- expand.grid(x = 1:8, y = 1:8)
  mu <- exp(1.5 + 0.8 * as.numeric(scale(grid$x)))
  counts <- rbind(stats::rpois(64, mu) + 1, stats::rpois(64, mu) + 1)
  spe <- grid_spe(counts, side = 8)
  fit <- function(two_part) {
    SummarizedExperiment::rowData(tweedspot(
      spe, two_part = two_part, smooth_k = 10, size_factor_method = "library_size",
      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
    ))$tweedspot_pval
  }
  expect_equal(fit(TRUE), fit(FALSE))
})

test_that("flag mode keeps everything and matches drop mode", {
  set.seed(2)
  counts <- rbind(sparse_counts(3), matrix(stats::rpois(4 * 36, 5), 4))
  counts[, 36] <- 0
  spe <- grid_spe(counts)
  args <- list(input = spe, smooth_k = 10, exclude = "f5", outlier_nmads = 3)
  dropped <- do.call(filter_features, args)
  flagged <- do.call(filter_features, c(args, action = "flag"))

  expect_identical(dim(flagged), dim(spe))
  rd <- SummarizedExperiment::rowData(flagged)
  cd <- SummarizedExperiment::colData(flagged)
  expect_identical(rownames(spe)[rd$filter_keep], rownames(dropped))
  expect_identical(colnames(spe)[cd$filter_keep], colnames(dropped))
  expect_identical(rd$filter_reason[rd$filter_keep], rep("", sum(rd$filter_keep)))
  expect_match(rd$filter_reason[1], "low_model_support")
  expect_identical(rd$filter_reason[5], "excluded")
  expect_match(cd$filter_reason[36], "low_total|outlier|empty")
  expect_identical(S4Vectors::metadata(flagged)$filter_features$action, "flag")
})

test_that("multiple failed criteria are all listed", {
  spe <- grid_spe(rbind(sparse_counts(2, value = 1), rep(4, 36)))
  flagged <- filter_features(spe, smooth_k = 10, action = "flag")
  reason <- SummarizedExperiment::rowData(flagged)$filter_reason[1]
  expect_identical(strsplit(reason, ";")[[1]], c("low_total", "low_model_support"))
})

test_that("location QC thresholds remove low-quality locations", {
  counts <- matrix(rep(4:9, length.out = 3 * 36), nrow = 3, ncol = 36)
  counts[, 1] <- c(1, 0, 0)        # low total, one feature detected
  counts[, 2] <- c(2, 2, 0)        # two features detected
  spe <- grid_spe(counts)
  base <- list(input = spe, min_total = 0, min_pct_locations = 0, model_aware = FALSE)
  kept <- function(...) colnames(do.call(filter_features, c(base, list(...))))

  expect_false("loc1" %in% kept(min_total_per_location = 2))
  expect_true("loc2" %in% kept(min_total_per_location = 2))
  expect_false(any(c("loc1", "loc2") %in% kept(min_features_per_location = 3)))
  out <- kept(outlier_nmads = 3)
  expect_false("loc1" %in% out)
  expect_identical(length(out), 35L)
  expect_length(kept(), 36L)
})

test_that("outlier detection is one-sided and robust to constant data", {
  x <- c(1, 9, 10, 11, 12, 30)
  expect_identical(tweedspot_low_outlier(x, 3),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  # A high value is never flagged, however extreme.
  expect_false(tweedspot_low_outlier(c(x, 1e6), 3)[7])
  expect_false(any(tweedspot_low_outlier(rep(5, 10), 3)))
})

test_that("max_excluded_fraction filters on signal from excluded features", {
  counts <- rbind(rep(10, 36), rep(10, 36), rep(1, 36))
  counts[3, 1] <- 30               # loc1: 30 of 50 counts from the excluded feature
  spe <- grid_spe(counts)
  rownames(spe) <- c("GENEA", "GENEB", "MT-CO1")
  res <- filter_features(spe, min_total = 0, min_pct_locations = 0, model_aware = FALSE,
                         exclude_mito = TRUE, max_excluded_fraction = 0.2)
  expect_false("loc1" %in% colnames(res))
  expect_identical(ncol(res), 35L)
  expect_false("MT-CO1" %in% rownames(res))
  # Works the same way for control probes selected by pattern.
  rownames(spe) <- c("GENEA", "GENEB", "NegPrb1")
  res2 <- filter_features(spe, min_total = 0, min_pct_locations = 0, model_aware = FALSE,
                          exclude_pattern = "^NegPrb", max_excluded_fraction = 0.2)
  expect_identical(colnames(res2), colnames(res))
  expect_error(filter_features(spe, max_excluded_fraction = 2), "max_excluded_fraction")
})

test_that("location_mask accepts logical or 0/1 columns and validates input", {
  spe <- grid_spe(matrix(5, nrow = 2, ncol = 36))
  SummarizedExperiment::colData(spe)$in_tissue <- rep(c(1, 0), 18)
  res <- filter_features(spe, min_total = 0, min_pct_locations = 0,
                         model_aware = FALSE, location_mask = "in_tissue")
  expect_identical(ncol(res), 18L)
  expect_error(filter_features(spe, location_mask = "missing"), "location_mask")
  SummarizedExperiment::colData(spe)$label <- rep("a", 36)
  expect_error(filter_features(spe, location_mask = "label"), "logical")
})
