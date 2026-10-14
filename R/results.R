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
#' @param group For results of `tweedspot(group = ...)`: `NULL` (default)
#'   returns the usual per-feature `rowData()` results; `TRUE` returns the full
#'   feature-by-group table; one or more group labels return only those
#'   groups' rows.
#'
#' @return A `data.frame`. Without `group`, one row per feature with columns
#'   `feature` (row name), `feature_name` (when `rowData()` has a
#'   `feature_name` or `gene_name` column that differs from the row names),
#'   `stat` (test statistic of the spatial smooth), `pval`, `padj`, `edf`
#'   (effective degrees of freedom of the smooth), `dev_expl` (deviance
#'   explained), and, when available, `p` (Tweedie variance power) and
#'   `family`. With `group`, one row per feature and group with columns `feature`, `group`, `n_locations`, `n_detected`,
#'   `tested`, `reason` (why a pair was not tested), `stat`, `pval`, `padj`,
#'   `padj_feature` (stage-wise screening p-value), `edf`, `dev_expl`, `p`, and
#'   `family`.
#'
#' @examples
#' counts <- matrix(c(1, 2, 3, 4, 2, 2), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("gene1", "gene2")
#' colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = counts),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(counts),
#'     gene_name = rownames(counts),
#'     tweedspot_stat = c(4, 2),
#'     tweedspot_pval = c(0.01, 0.05),
#'     tweedspot_padj = c(0.02, 0.05),
#'     tweedspot_edf = c(2, 2),
#'     tweedspot_dev_expl = c(0.3, 0.1)
#'   ),
#'   colData = S4Vectors::DataFrame(sample_id = rep("sample1", ncol(counts))),
#'   spatialCoords = cbind(x = seq_len(ncol(counts)), y = seq_len(ncol(counts)))
#' )
#' get_tweedspot_results(spe)
#' @export
get_tweedspot_results <- function(input,
                                  sort_by = c("padj", "none", "pval", "stat", "dev_expl"),
                                  decreasing = NULL,
                                  group = NULL) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  sort_by <- match.arg(sort_by)
  if (!is.null(group)) {
    return(tweedspot_group_results(input, group, sort_by, decreasing))
  }
  tweedspot_check_locations(input, S4Vectors::metadata(input)$tweedspot$locations)
  rd <- SummarizedExperiment::rowData(input)
  required <- paste0("tweedspot_", c("stat", "pval", "padj", "edf", "dev_expl"))
  missing_cols <- setdiff(required, colnames(rd))
  if (length(missing_cols)) {
    stop(sprintf(
      "TweedSpot result columns not found in `rowData(input)`: %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # A clean results table: one row per feature, short column names.
  feature_ids <- rownames(input)
  if (is.null(feature_ids)) {
    feature_ids <- paste0("feature", seq_len(nrow(input)))
  }
  res <- data.frame(feature = feature_ids, stringsAsFactors = FALSE)
  name_col <- tweedspot_name_cols(input)
  if (length(name_col)) {
    names_vec <- as.character(rd[[name_col[1L]]])
    if (!identical(names_vec, feature_ids)) {
      res$feature_name <- names_vec
    }
  }
  result_cols <- c("stat", "pval", "padj", "edf", "dev_expl", "p", "family")
  for (col in result_cols) {
    src <- paste0("tweedspot_", col)
    if (src %in% colnames(rd)) {
      res[[col]] <- rd[[src]]
    }
  }
  if (sort_by == "none") {
    return(res)
  }

  if (is.null(decreasing)) {
    decreasing <- sort_by %in% c("stat", "dev_expl")
  }
  if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
    stop("`decreasing` must be `TRUE`, `FALSE`, or `NULL`.")
  }
  res <- tweedspot_order_results(res, key_col = sort_by, decreasing = decreasing)
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
#' @param group For results of `tweedspot(group = ...)`: a group label (or
#'   several) to return the top hits within those groups, or `TRUE` for the top
#'   feature-group pairs overall. Default `NULL` uses the per-feature results.
#'
#' @return A ranked `data.frame` containing the top TweedSpot hits.
#'
#' @examples
#' counts <- matrix(c(1, 2, 3, 4, 2, 2), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("gene1", "gene2")
#' colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = counts),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(counts),
#'     gene_name = rownames(counts),
#'     tweedspot_stat = c(4, 2),
#'     tweedspot_pval = c(0.01, 0.05),
#'     tweedspot_padj = c(0.02, 0.05),
#'     tweedspot_edf = c(2, 2),
#'     tweedspot_dev_expl = c(0.3, 0.1)
#'   ),
#'   colData = S4Vectors::DataFrame(sample_id = rep("sample1", ncol(counts))),
#'   spatialCoords = cbind(x = seq_len(ncol(counts)), y = seq_len(ncol(counts)))
#' )
#' top_spatial_genes(spe, n = 1)
#' @export
top_spatial_genes <- function(input,
                              n = 10L,
                              rank_by = c("padj", "pval", "stat", "dev_expl"),
                              group = NULL) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  rank_by <- match.arg(rank_by)
  tweedspot_check_scalar(n, "n", min_value = 1)
  n <- as.integer(n)
  res <- get_tweedspot_results(input, sort_by = rank_by, group = group)
  if (!is.null(group)) {
    res <- res[res$tested, , drop = FALSE]
    if (!isTRUE(group)) {
      untested <- setdiff(as.character(group), unique(res$group))
      if (length(untested)) {
        message("No features were tested in group(s): ", paste(untested, collapse = ", "),
                " (see the `reason` column of get_tweedspot_results(group = TRUE)).")
      }
    }
    if (!isTRUE(group)) {
      # Top n within each requested group.
      res <- do.call(rbind, lapply(split(res, factor(res$group, levels = as.character(group))),
                                   function(d) d[seq_len(min(n, nrow(d))), , drop = FALSE]))
      rownames(res) <- NULL
      return(res)
    }
  }
  res[seq_len(min(n, nrow(res))), , drop = FALSE]
}

tweedspot_group_results <- function(input, group, sort_by, decreasing) {
  stored <- S4Vectors::metadata(input)$tweedspot_group
  if (is.null(stored)) {
    stop("No group-specific TweedSpot results found; run `tweedspot(group = ...)` first.")
  }
  res <- tweedspot_align_group_results(input, stored)
  if (!isTRUE(group)) {
    if (!is.character(group) && !is.factor(group) && !is.numeric(group) && !is.logical(group)) {
      stop("`group` must be NULL, TRUE, or one or more group labels.")
    }
    group <- as.character(group)
    unknown <- setdiff(group, unique(res$group))
    if (length(unknown)) {
      stop("Unknown group label(s): ", paste(unknown, collapse = ", "),
           ". Available: ", paste(unique(res$group), collapse = ", "), ".")
    }
    res <- res[res$group %in% group, , drop = FALSE]
  }
  if (sort_by != "none") {
    if (is.null(decreasing)) {
      decreasing <- sort_by %in% c("stat", "dev_expl")
    }
    if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
      stop("`decreasing` must be `TRUE`, `FALSE`, or `NULL`.")
    }
    key <- res[[sort_by]]
    stat <- res$stat
    ord <- if (decreasing) {
      order(is.na(key), -key, is.na(stat), -stat)
    } else {
      order(is.na(key), key, is.na(stat), -stat)
    }
    res <- res[ord, , drop = FALSE]
  }
  rownames(res) <- NULL
  res
}

# Keep the group table in step with the current object: drop rows for features
# that were removed by subsetting, re-index the rest, and warn if locations
# changed after the run.
tweedspot_align_group_results <- function(input, stored) {
  res <- stored$results
  tweedspot_check_locations(input, stored$locations)
  current <- rownames(input)
  if (identical(current, stored$features)) {
    return(res)
  }
  if (!is.null(current) && !anyDuplicated(current) && !is.null(stored$features) &&
      !anyDuplicated(stored$features) && all(current %in% stored$features)) {
    res <- res[res$feature %in% current, , drop = FALSE]
    res$feature_index <- match(res$feature, current)
    return(res)
  }
  warning("Features changed after `tweedspot()` in a way the stored group results ",
          "cannot be matched to; re-run `tweedspot()` on the current object.", call. = FALSE)
  res
}
