#' Filter low-information genes before spatial testing
#'
#' `filter_genes()` removes genes with too few total counts or too little
#' spatial prevalence to support stable downstream modeling.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a count assay.
#' @param assay_name Name of the assay to filter. Default `"counts"`.
#' @param filter_genes_ncounts Minimum total count across all spatial locations
#'   required to retain a gene. Default `5`.
#' @param filter_genes_pcspots Minimum percent of spatial locations with
#'   non-zero counts required to retain a gene. Default `5`.
#' @param filter_genes_nspots Minimum number of spatial locations with non-zero
#'   counts required to retain a gene. Default `0`.
#' @param filter_genes_mean Minimum mean count across spatial locations required
#'   to retain a gene. Default `0`.
#' @param filter_genes_var Minimum variance across spatial locations required to
#'   retain a gene. Default `0`.
#' @param exclude_mito Logical; if `TRUE`, remove mitochondrial genes using
#'   `rowData(input)` annotations or gene-name heuristics. Default `FALSE`.
#' @param exclude_ribo Logical; if `TRUE`, remove ribosomal genes using
#'   `rowData(input)` annotations or gene-name heuristics. Default `FALSE`.
#' @param gene_name_col Optional `rowData(input)` column containing gene names
#'   to use for mitochondrial/ribosomal exclusion. If `NULL`, row names are
#'   used. Default `NULL`.
#' @param gene_biotype_col Optional `rowData(input)` column containing gene
#'   biotypes or categories to use for mitochondrial/ribosomal exclusion.
#'   Default `NULL`.
#' @param drop_zero_spots Logical; if `TRUE`, drop spatial locations with zero
#'   total counts after gene filtering. Default `TRUE`.
#' @param report Logical; if `TRUE`, print a short summary of how many genes and
#'   spots were retained or removed. Default `FALSE`.
#'
#' @return A filtered `SpatialExperiment` containing only genes that pass all
#'   requested thresholds and, by default, only spatial locations with positive
#'   total counts after filtering.
#'
#' @examples
#' \dontrun{
#' library(STexampleData)
#' spe <- ST_mouseOB()
#' spe <- filter_genes(
#'   spe,
#'   filter_genes_ncounts = 10,
#'   filter_genes_pcspots = 5,
#'   filter_genes_nspots = 3,
#'   drop_zero_spots = TRUE
#' )
#' }
#' @export
filter_genes <- function(input,
                         assay_name = "counts",
                         filter_genes_ncounts = 5,
                         filter_genes_pcspots = 5,
                         filter_genes_nspots = 0,
                         filter_genes_mean = 0,
                         filter_genes_var = 0,
                         exclude_mito = FALSE,
                         exclude_ribo = FALSE,
                         gene_name_col = NULL,
                         gene_biotype_col = NULL,
                         drop_zero_spots = TRUE,
                         report = FALSE) {

  stopifnot(methods::is(input, "SpatialExperiment"))
  tweedspot_check_scalar(filter_genes_ncounts, "filter_genes_ncounts", min_value = 0)
  tweedspot_check_scalar(filter_genes_pcspots, "filter_genes_pcspots", min_value = 0,
                         max_value = 100)
  tweedspot_check_scalar(filter_genes_nspots, "filter_genes_nspots", min_value = 0)
  tweedspot_check_scalar(filter_genes_mean, "filter_genes_mean", min_value = 0)
  tweedspot_check_scalar(filter_genes_var, "filter_genes_var", min_value = 0)
  if (!is.logical(exclude_mito) || length(exclude_mito) != 1L || is.na(exclude_mito)) {
    stop("`exclude_mito` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(exclude_ribo) || length(exclude_ribo) != 1L || is.na(exclude_ribo)) {
    stop("`exclude_ribo` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(drop_zero_spots) || length(drop_zero_spots) != 1L || is.na(drop_zero_spots)) {
    stop("`drop_zero_spots` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("`report` must be `TRUE` or `FALSE`.")
  }

  Y <- as.matrix(SummarizedExperiment::assay(input, assay_name))
  keep_ncounts <- rowSums(Y) >= filter_genes_ncounts
  keep_pcspots <- rowMeans(Y > 0) * 100 >= filter_genes_pcspots
  keep_nspots <- rowSums(Y > 0) >= filter_genes_nspots
  keep_mean <- rowMeans(Y) >= filter_genes_mean
  keep_var <- apply(Y, 1, stats::var) >= filter_genes_var
  keep <- keep_ncounts & keep_pcspots & keep_nspots & keep_mean & keep_var
  n_genes_before <- nrow(input)
  n_spots_before <- ncol(input)

  if (exclude_mito || exclude_ribo) {
    drop_genes <- tweedspot_annotation_exclusions(
      input = input,
      exclude_mito = exclude_mito,
      exclude_ribo = exclude_ribo,
      gene_name_col = gene_name_col,
      gene_biotype_col = gene_biotype_col
    )
    keep <- keep & !drop_genes
  }

  output <- input[keep, ]
  dropped_zero_spots <- 0L
  if (drop_zero_spots) {
    output_counts <- as.matrix(SummarizedExperiment::assay(output, assay_name))
    keep_spots <- colSums(output_counts) > 0
    dropped_zero_spots <- sum(!keep_spots)
    output <- output[, keep_spots]
  }
  if (report) {
    criteria_removed <- c(
      ncounts = sum(!keep_ncounts),
      pcspots = sum(!keep_pcspots),
      nspots = sum(!keep_nspots),
      mean = sum(!keep_mean),
      variance = sum(!keep_var)
    )
    if (exclude_mito || exclude_ribo) {
      criteria_removed <- c(criteria_removed, annotations = sum(drop_genes))
    }
    message(
      "filter_genes retained ", nrow(output), "/", n_genes_before, " genes and ",
      ncol(output), "/", n_spots_before, " spots."
    )
    message(
      "Genes removed by criterion: ",
      paste(names(criteria_removed), criteria_removed, sep = "=", collapse = ", "),
      ". Spots dropped due to zero library size after filtering: ",
      dropped_zero_spots, "."
    )
  }
  output
}
