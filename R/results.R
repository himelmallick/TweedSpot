#' Extract TweedSpot results as a tidy data frame
#'
#' `get_tweedspot_results()` pulls the TweedSpot result columns from
#' `rowData(input)` into a base `data.frame` and can optionally sort the output.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] containing TweedSpot
#'   results in `rowData(input)`.
#' @param sort_by Column to sort by. Options are `"none"`, `"padj"`, `"pval"`,
#'   `"stat"`, or `"dev_expl"`. Default `"padj"`.
#' @param decreasing Logical; override the default sort direction for
#'   `sort_by`. If `NULL`, p-values sort ascending and statistics/effect-size
#'   columns sort descending.
#'
#' @return A `data.frame` containing the original `rowData(input)` columns and
#'   TweedSpot result columns.
#'
#' @examples
#' \dontrun{
#' res <- get_tweedspot_results(spe_res)
#' head(res)
#' }
#' @export
get_tweedspot_results <- function(input,
                                  sort_by = c("padj", "none", "pval", "stat", "dev_expl"),
                                  decreasing = NULL) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  sort_by <- match.arg(sort_by)
  res <- as.data.frame(SummarizedExperiment::rowData(input))
  required <- c(
    "tweedspot_stat",
    "tweedspot_pval",
    "tweedspot_padj",
    "tweedspot_edf",
    "tweedspot_dev_expl"
  )
  missing_cols <- setdiff(required, colnames(res))
  if (length(missing_cols)) {
    stop(sprintf(
      "TweedSpot result columns not found in `rowData(input)`: %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (sort_by == "none") {
    return(res)
  }

  default_decreasing <- switch(
    sort_by,
    padj = FALSE,
    pval = FALSE,
    stat = TRUE,
    dev_expl = TRUE
  )
  if (is.null(decreasing)) {
    decreasing <- default_decreasing
  }
  if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
    stop("`decreasing` must be `TRUE`, `FALSE`, or `NULL`.")
  }

  key_col <- switch(
    sort_by,
    padj = "tweedspot_padj",
    pval = "tweedspot_pval",
    stat = "tweedspot_stat",
    dev_expl = "tweedspot_dev_expl"
  )
  res <- tweedspot_order_results(res, key_col = key_col, decreasing = decreasing)
  rownames(res) <- NULL
  res
}

#' Return the top TweedSpot hits
#'
#' `top_spatial_genes()` returns the top `n` rows from
#' [get_tweedspot_results()] using a user-selected ranking metric.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] containing TweedSpot
#'   results in `rowData(input)`.
#' @param n Number of rows to return. Default `10`.
#' @param rank_by Ranking column. Options are `"padj"`, `"pval"`, `"stat"`, or
#'   `"dev_expl"`. Default `"padj"`.
#'
#' @return A ranked `data.frame` containing the top TweedSpot hits.
#'
#' @examples
#' \dontrun{
#' top_spatial_genes(spe_res, n = 10)
#' }
#' @export
top_spatial_genes <- function(input,
                              n = 10L,
                              rank_by = c("padj", "pval", "stat", "dev_expl")) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  rank_by <- match.arg(rank_by)
  tweedspot_check_scalar(n, "n", min_value = 1)
  n <- as.integer(n)
  res <- get_tweedspot_results(input, sort_by = rank_by)
  res[seq_len(min(n, nrow(res))), , drop = FALSE]
}
