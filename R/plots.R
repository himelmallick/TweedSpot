#' Plot observed spatial expression for a gene
#'
#' `plot_spatial_gene()` visualizes observed expression for a single gene across
#' spatial coordinates.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a count assay and
#'   spatial coordinates.
#' @param gene Gene identifier, row name, or row index to plot.
#' @param assay_name Name of the assay to visualize. Default `"counts"`.
#' @param transform Transformation to apply to observed counts. Options are
#'   `"log1p"`, `"sqrt"`, or `"none"`. Default `"log1p"`.
#' @param point_size Point size for the scatter plot. Default `1.8`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' counts <- matrix(c(0, 1, 2, 3, 1, 1, 1, 1), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("gene1", "gene2")
#' colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = counts),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(counts),
#'     gene_name = rownames(counts)
#'   ),
#'   colData = S4Vectors::DataFrame(batch = rep(c("A", "B"), length.out = ncol(counts))),
#'   spatialCoords = cbind(x = seq_len(ncol(counts)), y = seq_len(ncol(counts)))
#' )
#' plot_spatial_gene(spe, gene = "gene1")
#' @export
plot_spatial_gene <- function(input,
                              gene,
                              assay_name = "counts",
                              transform = c("log1p", "sqrt", "none"),
                              point_size = 1.8) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_spatial_gene()` requires the `ggplot2` package.")
  }
  transform <- match.arg(transform)
  tweedspot_check_scalar(point_size, "point_size", min_value = 0)

  gene_info <- tweedspot_resolve_gene(input, gene)
  coords <- tweedspot_plot_coords(input)
  expr <- as.numeric(SummarizedExperiment::assay(input, assay_name)[gene_info$index, ])
  value <- switch(
    transform,
    log1p = log1p(expr),
    sqrt = sqrt(expr),
    none = expr
  )
  legend_title <- switch(
    transform,
    log1p = "log1p(counts)",
    sqrt = "sqrt(counts)",
    none = "counts"
  )

  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    value = value
  )

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = value)) +
    ggplot2::geom_point(size = point_size, alpha = 0.95) +
    ggplot2::scale_color_gradientn(
      colors = grDevices::hcl.colors(9, "YlOrRd", rev = TRUE),
      name = legend_title
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    tweedspot_plot_theme() +
    ggplot2::labs(
      title = paste0("Observed spatial expression: ", gene_info$label),
      x = "x",
      y = "y"
    )
}

#' Plot the fitted spatial signal for a gene
#'
#' `plot_spatial_fit()` refits the TweedSpot gene-level GAM for a single gene
#' and visualizes the fitted spatial smooth contribution across spatial
#' locations.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a count assay and
#'   spatial coordinates.
#' @param gene Gene identifier, row name, or row index to plot.
#' @param assay_name Name of the assay to model. Default `"counts"`.
#' @param covariates Optional one-sided formula specifying columns from
#'   `colData(input)` to adjust for, for example `~ batch + age + sex`.
#'   Character vectors of `colData(input)` column names are also supported.
#' @param family `mgcv` family for the count model. Default `"tw"`.
#' @param fit_method Smoothing parameter estimation method passed to `mgcv`.
#'   Default `"REML"`.
#' @param use_bam Logical; use [mgcv::bam()] instead of [mgcv::gam()]. Default
#'   `FALSE`.
#' @param bam_discrete Logical; when `use_bam = TRUE`, use discretized fitting.
#'   Default `TRUE`.
#' @param bam_nthreads Integer; number of threads for `mgcv::bam()`. Default 1.
#' @param smooth_k Optional basis dimension for the spatial smooth. Default
#'   `NULL`.
#' @param point_size Point size for the scatter plot. Default `1.8`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' counts <- matrix(
#'   c(0, 1, 2, 3, 4, 5, 5, 4, 5, 4, 3, 2),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("gene1", "gene2")
#' colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = counts),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(counts),
#'     gene_name = rownames(counts)
#'   ),
#'   colData = S4Vectors::DataFrame(batch = rep(c("A", "B"), length.out = ncol(counts))),
#'   spatialCoords = cbind(x = seq_len(ncol(counts)), y = seq_len(ncol(counts)))
#' )
#' plot_spatial_fit(spe, gene = "gene1", covariates = ~ batch, family = "poisson")
#' @export
plot_spatial_fit <- function(input,
                             gene,
                             assay_name = "counts",
                             covariates = NULL,
                             family = "tw",
                             fit_method = "REML",
                             use_bam = FALSE,
                             bam_discrete = TRUE,
                             bam_nthreads = 1L,
                             smooth_k = NULL,
                             point_size = 1.8) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_spatial_fit()` requires the `ggplot2` package.")
  }
  tweedspot_check_scalar(point_size, "point_size", min_value = 0)

  fit_data <- tweedspot_prepare_plot_fit_data(
    input = input,
    gene = gene,
    assay_name = assay_name,
    covariates = covariates,
    family = family,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    bam_nthreads = bam_nthreads,
    smooth_k = smooth_k
  )

  ggplot2::ggplot(fit_data$df, ggplot2::aes(x = x, y = y, color = value)) +
    ggplot2::geom_point(size = point_size, alpha = 0.95) +
    ggplot2::scale_color_gradient2(
      low = "#1B4D6B",
      mid = "#F6F5F2",
      high = "#A33D2B",
      midpoint = 0,
      name = "fitted smooth"
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    tweedspot_plot_theme() +
    ggplot2::labs(
      title = paste0("Fitted spatial signal: ", fit_data$label),
      x = "x",
      y = "y"
    )
}
