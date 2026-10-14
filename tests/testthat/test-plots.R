plot_spe <- function(fit = TRUE) {
  set.seed(5)
  grid <- expand.grid(x = 1:10, y = 1:10)
  counts <- rbind(
    band = stats::rpois(100, ifelse(grid$y > 5, 12, 2)),
    blob = stats::rpois(100, 2 + 10 * exp(-((grid$x - 3)^2 + (grid$y - 3)^2) / 5)),
    flat = stats::rpois(100, 5)
  )
  colnames(counts) <- paste0("loc", seq_len(ncol(counts)))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(feature_name = c("Band", "Blob", "Flat")),
    spatialCoords = as.matrix(grid)
  )
  if (fit) {
    spe <- tweedspot(spe, smooth_k = 8, size_factor_method = "library_size",
                     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  }
  spe
}

test_that("tweedspot() stores the settings that plotting reuses", {
  spe <- plot_spe()
  stored <- S4Vectors::metadata(spe)$tweedspot
  expect_identical(stored$smooth_k, 8)
  expect_identical(stored$size_factor_method, "library_size")
  settings <- tweedspot_plot_settings(spe)
  expect_identical(settings$smooth_k, 8)
  expect_identical(tweedspot_plot_settings(spe, list(smooth_k = 5))$smooth_k, 5)
  expect_error(tweedspot_plot_settings(spe, list(bogus = 1)), "Unknown model setting")
})

test_that("plot_feature draws three panels from the stored model", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spe <- plot_spe()
  p <- plot_feature(spe, "band")
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 2L)
  expect_identical(p$patches$annotation$title, "Band (band)")
  expect_s3_class(plot_feature(spe, "Blob", transform = "none", flip_y = FALSE,
                               point_size = 1), "patchwork")
  expect_s3_class(plot_feature(spe, 3, smooth_k = 5), "patchwork")
})

test_that("plotted fits match the tested model", {
  spe <- plot_spe()
  fit <- tweedspot_plot_fits(spe, 1L, tweedspot_plot_settings(spe))[[1]]
  expect_length(fit$effect, ncol(spe))
  expect_true(all(fit$fitted > 0))
  # The smooth is centred, so the log2 effect averages to about zero.
  expect_lt(abs(mean(fit$effect)), 1e-6)
  # The band feature is higher where y > 5.
  y <- SpatialExperiment::spatialCoords(spe)[, 2]
  expect_gt(mean(fit$effect[y > 5]), mean(fit$effect[y <= 5]))
})

test_that("plot_feature validates its inputs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spe <- plot_spe()
  expect_error(plot_feature(spe, 99), "out of bounds")
  expect_error(plot_feature(spe, "missing"), "was not found")
  expect_error(plot_feature(spe, "band", point_size = -1), "point_size")
  expect_error(plot_feature(spe, "band", formats = "jpg"), "formats")
  expect_error(plot_feature(spe, "band", output_dir = c("a", "b")), "output_dir")
  expect_error(plot_feature(spe, "band", smooth_kk = 5), "Unknown model setting")
})

test_that("plot_feature saves to a new output directory", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spe <- plot_spe()
  dir <- file.path(tempfile("tweedspot-plots"), "nested")
  p <- plot_feature(spe, "band", output_dir = dir, formats = c("pdf", "png"), dpi = 72)
  files <- attr(p, "files")
  expect_true(dir.exists(dir))
  expect_identical(basename(files), c("band_tweedspot.pdf", "band_tweedspot.png"))
  expect_true(all(file.exists(files) & file.size(files) > 0))
  unlink(dirname(dir), recursive = TRUE)
})

test_that("plot_top_features honours n_top, features, labels and title", {
  skip_if_not_installed("ggplot2")
  spe <- plot_spe()
  panels <- function(p) levels(p$data$panel)

  p2 <- plot_top_features(spe, n_top = 2)
  expect_s3_class(p2, "ggplot")
  expect_length(panels(p2), 2L)
  expect_match(panels(p2)[1], "^(Band \\(band\\)|Blob \\(blob\\))\\nadj\\. p (<|=) .*% deviance$")

  p_sel <- plot_top_features(spe, features = c("flat", "Band"), label_stats = FALSE)
  expect_identical(panels(p_sel), c("Flat (flat)", "Band (band)"))

  expect_null(p2$labels$title)
  expect_identical(plot_top_features(spe, n_top = 1, title = "SVFs")$labels$title, "SVFs")
  expect_message(p_all <- plot_top_features(spe, n_top = 10), "Only 3")
  expect_length(panels(p_all), 3L)
})

test_that("plot_top_features ranks by adjusted p-value", {
  spe <- plot_spe()
  rd <- SummarizedExperiment::rowData(spe)
  rd$tweedspot_padj <- c(0.5, 0.01, 0.2)
  SummarizedExperiment::rowData(spe) <- rd
  expect_identical(tweedspot_top_features_ordinary(spe, 3, NULL)$index, c(2L, 3L, 1L))
  rd$tweedspot_padj[2] <- NA
  SummarizedExperiment::rowData(spe) <- rd
  expect_identical(tweedspot_top_features_ordinary(spe, 3, NULL)$index, c(3L, 1L))
})

test_that("plot_top_features validates inputs and saves files", {
  skip_if_not_installed("ggplot2")
  spe <- plot_spe()
  expect_error(plot_top_features(plot_spe(fit = FALSE)), "run `tweedspot\\(\\)` first")
  expect_error(plot_top_features(spe, n_top = 1.5), "whole number")
  expect_error(plot_top_features(spe, n_top = 0), "n_top")
  expect_error(plot_top_features(spe, title = 1), "title")
  dir <- tempfile("tweedspot-top")
  p <- plot_top_features(spe, n_top = 2, output_dir = dir, formats = "png", dpi = 72)
  expect_identical(basename(attr(p, "files")), "tweedspot_top_features.png")
  expect_true(file.exists(attr(p, "files")))
  unlink(dir, recursive = TRUE)
})

test_that("plot helpers format values and size markers sensibly", {
  expect_identical(tweedspot_format_p(0), "< 1e-300")
  expect_identical(tweedspot_format_p(0.0123), "= 0.012")
  expect_identical(tweedspot_format_p(1.5e-12), "= 1.5e-12")
  expect_identical(tweedspot_format_p(NA_real_), "NA")
  expect_identical(tweedspot_file_stem("MT-CO1 / x"), "MT-CO1_x")
  expect_identical(tweedspot_file_stem("///"), "feature")
  coords <- as.matrix(expand.grid(x = 1:20, y = 1:20))
  small <- tweedspot_auto_point_size(coords, 1, 1)
  large <- tweedspot_auto_point_size(coords, 4, 4)
  expect_true(is.finite(small) && small > 0)
  expect_gt(large, small)
  expect_equal(tweedspot_auto_point_size(coords[1, , drop = FALSE], 1, 1), 2)
  expect_identical(scales_squish(c(-2, 0.5, 3), c(0, 1)), c(0, 0.5, 1))
})

test_that("plot colours come from the Mako and Rocket colour maps", {
  mako <- toupper(grDevices::hcl.colors(9, "Mako", rev = TRUE))
  rocket <- toupper(grDevices::hcl.colors(9, "Rocket", rev = TRUE))
  expect_identical(toupper(tweedspot_seq_colours()), mako)
  div <- toupper(tweedspot_div_colours())
  expect_length(div, 11L)
  expect_identical(div[1:5], rev(mako[2:6]))
  expect_identical(div[7:11], rocket[2:6])
})
