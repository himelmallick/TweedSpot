#' Filter low-information features and low-quality locations
#'
#' `filter_features()` prepares a spatial dataset for [tweedspot()]. It removes
#' (or flags) features with too little signal to model, features that cannot
#' support the spatial smooth TweedSpot will fit, named features such as
#' control probes, and low-quality spatial locations. All criteria are
#' modality agnostic and apply to genes, proteins, peaks, or metabolites, and
#' to spots, cells, or pixels; mitochondrial and ribosomal gene removal is
#' available as a transcriptomics-specific option.
#'
#' Filtering happens in this order:
#'
#' 1. **Named exclusions** (`exclude`, `exclude_pattern`, `exclude_mito`,
#'    `exclude_ribo`) mark features to remove.
#' 2. **Location QC** is computed on the remaining features: minimum total
#'    signal, minimum detected features, lower-tail outliers, the fraction of
#'    signal from excluded features, and an optional mask.
#' 3. **Feature thresholds** are computed on the locations that pass QC,
#'    including the model-aware minimum below.
#' 4. Locations with no remaining signal are removed (`drop_empty_locations`).
#'
#' **Model-aware minimum.** With `model_aware = TRUE`, a feature must be
#' detected in at least as many locations as the spatial smooth has basis
#' functions (`smooth_k`, or TweedSpot's default of up to 30). Features below
#' this cannot support the smooth and would otherwise fail or give unstable
#' fits. Pass the same `smooth_k` you will give [tweedspot()]. With `group`,
#' the minimum is checked within each group using that group's locations, and
#' a feature is kept if at least one group can support it; the per-group
#' result is stored in `metadata(output)$filter_features$group_support`. For the
#' two-part model the presence/absence part additionally needs that many
#' detected *and* undetected locations; features that fail only this are kept,
#' and [tweedspot()] then tests them with the count part alone.
#'
#' All feature criteria depend only on abundance and detection, not on spatial
#' arrangement, so filtering does not bias the spatial test's p-values.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a non-negative
#'   feature-by-location assay.
#' @param assay_name Name of the assay to filter. Default `"counts"`.
#' @param min_total Minimum total signal (sum across locations) for a feature.
#'   Default `5`.
#' @param min_pct_locations Minimum percent of locations in which the feature
#'   is detected. Default `5`.
#' @param min_locations Minimum number of locations in which the feature is
#'   detected. Default `0`.
#' @param min_mean Minimum mean signal across locations. Default `0`.
#' @param min_var Minimum variance across locations. Default `0`.
#' @param detection_threshold A feature counts as detected at a location when
#'   its value is strictly greater than this. Default `0`, which suits counts;
#'   raise it for intensity data with a background floor.
#' @param model_aware Logical; require each feature to be detected in at least
#'   as many locations as the spatial smooth has basis functions (see
#'   Details). Default `TRUE`.
#' @param smooth_k Basis dimension of the spatial smooth, as passed to
#'   [tweedspot()]. Default `NULL` uses TweedSpot's default.
#' @param group Optional name of a discrete `colData(input)` column, as passed
#'   to `tweedspot(group = ...)`. The model-aware minimum is then checked within
#'   each group (of at least `min_group_size` locations), and a feature is
#'   removed only if no group can support it. Default `NULL`.
#' @param min_group_size Groups with fewer locations than this are ignored when
#'   `group` is given. Default `50`.
#' @param two_part Logical; whether you will fit [tweedspot()] with
#'   `two_part = TRUE`. Only affects reporting of features whose
#'   presence/absence part will be skipped. Default `FALSE`.
#' @param exclude Optional character vector of feature names to remove, matched
#'   exactly against row names and, if supplied, `feature_name_col`.
#' @param exclude_pattern Optional character vector of regular expressions;
#'   features whose row name or `feature_name_col` value matches any of them are
#'   removed, for example `"^NegPrb"`.
#' @param feature_name_col Optional `rowData(input)` column holding readable
#'   feature names, used by the exclusion arguments.
#' @param exclude_mito Logical; transcriptomics only. Remove mitochondrial
#'   genes (names beginning `MT-`/`mt-`, or a `feature_type_col` value
#'   containing "mitochond"). Default `FALSE`.
#' @param exclude_ribo Logical; transcriptomics only. Remove ribosomal protein
#'   genes (`RPL10A`, `RPLP0`, `RPSA`, `MRPL12`, ..., but not the `RPS6K`
#'   kinases), or a `feature_type_col` value containing "ribosom". Default
#'   `FALSE`.
#' @param feature_type_col Optional `rowData(input)` column of feature
#'   categories (for example gene biotypes) used by `exclude_mito` and
#'   `exclude_ribo`.
#' @param min_total_per_location Minimum total signal for a location. Default
#'   `0`.
#' @param min_features_per_location Minimum number of features detected at a
#'   location. Default `0`.
#' @param outlier_nmads Optional number of median absolute deviations; a
#'   location is removed when its log total signal or log number of detected
#'   features falls more than this many MADs below the median. Default `NULL`
#'   (off); `3` is a common choice.
#' @param max_excluded_fraction Optional maximum fraction of a location's
#'   signal that may come from excluded features (for example `0.2` with
#'   `exclude_mito = TRUE` for a mitochondrial-fraction filter, or with
#'   `exclude_pattern` for control-probe signal). Default `NULL` (off).
#' @param location_mask Optional name of a logical `colData(input)` column;
#'   locations where it is `FALSE` are removed, for example `"in_tissue"`.
#' @param drop_empty_locations Logical; remove locations with no signal after
#'   feature filtering. Default `TRUE`.
#' @param action `"drop"` (default) returns only retained features and
#'   locations. `"flag"` returns every feature and location, with
#'   `filter_keep` (logical) and `filter_reason` (failed criteria separated by
#'   `";"`, empty when kept) added to `rowData()` and `colData()`.
#' @param report Logical; print how many features and locations were kept and
#'   how many failed each criterion. Default `FALSE`.
#'
#' @return A `SpatialExperiment`. With `action = "drop"`, only retained
#'   features and locations; with `action = "flag"`, all of them with
#'   `filter_keep` and `filter_reason` columns. In both cases
#'   `S4Vectors::metadata(output)$filter_features` records the settings and the
#'   number of features and locations failing each criterion.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(stats::rpois(5 * 36, 4), nrow = 5)
#' counts[2, ] <- 0
#' counts[2, 1:2] <- 3
#' rownames(counts) <- c("GENE1", "RARE1", "MT-CO1", "NegPrb1", "GENE5")
#' colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'   assays = list(counts = counts),
#'   spatialCoords = as.matrix(expand.grid(x = 1:6, y = 1:6))
#' )
#' # Drop features, using a smooth with 10 basis functions
#' filter_features(spe, smooth_k = 10, exclude_pattern = "^NegPrb")
#' # Keep everything but record why each feature or location would be removed
#' flagged <- filter_features(spe, smooth_k = 10, exclude_mito = TRUE,
#'                            max_excluded_fraction = 0.3, action = "flag")
#' SummarizedExperiment::rowData(flagged)[, c("filter_keep", "filter_reason")]
#' @export
filter_features <- function(input,
                            assay_name = "counts",
                            min_total = 5,
                            min_pct_locations = 5,
                            min_locations = 0,
                            min_mean = 0,
                            min_var = 0,
                            detection_threshold = 0,
                            model_aware = TRUE,
                            smooth_k = NULL,
                            group = NULL,
                            min_group_size = 50,
                            two_part = FALSE,
                            exclude = NULL,
                            exclude_pattern = NULL,
                            feature_name_col = NULL,
                            exclude_mito = FALSE,
                            exclude_ribo = FALSE,
                            feature_type_col = NULL,
                            min_total_per_location = 0,
                            min_features_per_location = 0,
                            outlier_nmads = NULL,
                            max_excluded_fraction = NULL,
                            location_mask = NULL,
                            drop_empty_locations = TRUE,
                            action = c("drop", "flag"),
                            report = FALSE) {

  stopifnot(methods::is(input, "SpatialExperiment"))
  action <- match.arg(action)
  tweedspot_check_scalar(min_total, "min_total", min_value = 0)
  tweedspot_check_scalar(min_pct_locations, "min_pct_locations", min_value = 0,
                         max_value = 100)
  tweedspot_check_scalar(min_locations, "min_locations", min_value = 0)
  tweedspot_check_scalar(min_mean, "min_mean", min_value = 0)
  tweedspot_check_scalar(min_var, "min_var", min_value = 0)
  tweedspot_check_scalar(detection_threshold, "detection_threshold", min_value = 0)
  tweedspot_check_scalar(min_total_per_location, "min_total_per_location", min_value = 0)
  tweedspot_check_scalar(min_features_per_location, "min_features_per_location",
                         min_value = 0)
  if (!is.null(outlier_nmads)) {
    tweedspot_check_scalar(outlier_nmads, "outlier_nmads", min_value = 0)
  }
  if (!is.null(max_excluded_fraction)) {
    tweedspot_check_scalar(max_excluded_fraction, "max_excluded_fraction",
                           min_value = 0, max_value = 1)
  }
  if (!is.null(smooth_k)) {
    tweedspot_check_scalar(smooth_k, "smooth_k", min_value = 3)
  }
  tweedspot_check_scalar(min_group_size, "min_group_size", min_value = 1)
  labels <- if (is.null(group)) NULL else tweedspot_group_labels(input, group)
  for (flag in c("model_aware", "two_part", "exclude_mito", "exclude_ribo",
                 "drop_empty_locations", "report")) {
    tweedspot_check_flag(get(flag), flag)
  }

  Y <- tweedspot_get_assay(input, assay_name)
  if (min(Y) < 0) {
    stop("`filter_features()` expects a non-negative assay; `", assay_name,
         "` contains negative values.")
  }
  mask <- tweedspot_location_mask(input, location_mask)

  # 1. Named exclusions.
  excluded <- tweedspot_explicit_exclusions(
    input = input,
    exclude = exclude,
    exclude_pattern = exclude_pattern,
    feature_name_col = feature_name_col
  )
  mito_ribo <- rep(FALSE, nrow(input))
  if (exclude_mito || exclude_ribo) {
    mito_ribo <- tweedspot_annotation_exclusions(
      input = input,
      exclude_mito = exclude_mito,
      exclude_ribo = exclude_ribo,
      feature_name_col = feature_name_col,
      feature_type_col = feature_type_col
    )
  }
  named_out <- excluded | mito_ribo

  # 2. Location QC on the features that are not excluded by name.
  Y_measured <- Y[!named_out, , drop = FALSE]
  excluded_fraction <- NULL
  if (any(named_out)) {
    totals_all <- colSums(Y)
    excluded_fraction <- ifelse(
      totals_all > 0,
      colSums(Y[named_out, , drop = FALSE]) / totals_all,
      0
    )
  }
  location_fails <- tweedspot_location_filters(
    Y = Y_measured,
    detection_threshold = detection_threshold,
    min_total_per_location = min_total_per_location,
    min_features_per_location = min_features_per_location,
    outlier_nmads = outlier_nmads,
    excluded_fraction = excluded_fraction,
    max_excluded_fraction = max_excluded_fraction,
    mask = mask
  )
  keep_location <- !Reduce(`|`, location_fails)

  # 3. Feature thresholds on the locations that pass QC.
  coords <- SpatialExperiment::spatialCoords(input)[keep_location, , drop = FALSE]
  group_support <- NULL
  if (model_aware && !is.null(labels)) {
    # Support is checked within each group; the global minimum is then off.
    group_support <- tweedspot_group_support(
      Y[, keep_location, drop = FALSE], coords, labels[keep_location],
      smooth_k, min_group_size, detection_threshold
    )
    min_support <- 0L
  } else if (model_aware && sum(keep_location) > 0L) {
    min_support <- tweedspot_basis_k(coords, smooth_k)
  } else {
    min_support <- 0L
  }
  feature_fails <- tweedspot_feature_filters(
    Y = Y[, keep_location, drop = FALSE],
    min_total = min_total,
    min_pct_locations = min_pct_locations,
    min_locations = min_locations,
    min_mean = min_mean,
    min_var = min_var,
    detection_threshold = detection_threshold,
    min_support = min_support
  )
  if (!is.null(group_support)) {
    feature_fails$low_model_support <- if (ncol(group_support$supported)) {
      rowSums(group_support$supported) == 0L
    } else {
      rep(TRUE, nrow(input))
    }
  }
  feature_fails$excluded <- excluded
  feature_fails$mito_ribo <- mito_ribo
  keep_feature <- !Reduce(`|`, feature_fails)

  # 4. Locations with no remaining signal.
  location_fails$empty <- rep(FALSE, ncol(input))
  if (drop_empty_locations) {
    location_fails$empty <- keep_location &
      colSums(Y[keep_feature, , drop = FALSE]) <= 0
    keep_location <- keep_location & !location_fails$empty
  }

  binary_skipped <- 0L
  if (two_part && is.null(labels) && min_support > 0L && any(keep_feature)) {
    Yk <- Y[keep_feature, keep_location, drop = FALSE]
    n_det <- rowSums(Yk > 0)
    binary_skipped <- sum(n_det < min_support | (ncol(Yk) - n_det) < min_support)
  }

  feature_counts <- vapply(feature_fails, sum, integer(1))
  location_counts <- vapply(location_fails, sum, integer(1))
  summary_info <- list(
    action = action,
    n_features = c(before = nrow(input), kept = sum(keep_feature)),
    n_locations = c(before = ncol(input), kept = sum(keep_location)),
    feature_failures = feature_counts,
    location_failures = location_counts,
    min_support = if (is.null(group_support)) min_support else group_support$min_support,
    group_support = if (is.null(group_support)) NULL else group_support$supported,
    binary_part_skipped = binary_skipped,
    settings = list(
      assay_name = assay_name, min_total = min_total,
      min_pct_locations = min_pct_locations, min_locations = min_locations,
      min_mean = min_mean, min_var = min_var,
      detection_threshold = detection_threshold, model_aware = model_aware,
      smooth_k = smooth_k, group = group, min_group_size = min_group_size,
      two_part = two_part, exclude = exclude,
      exclude_pattern = exclude_pattern, exclude_mito = exclude_mito,
      exclude_ribo = exclude_ribo,
      min_total_per_location = min_total_per_location,
      min_features_per_location = min_features_per_location,
      outlier_nmads = outlier_nmads,
      max_excluded_fraction = max_excluded_fraction,
      location_mask = location_mask,
      drop_empty_locations = drop_empty_locations
    )
  )

  if (report) {
    tweedspot_report_filtering(
      n_features_kept = sum(keep_feature),
      n_features_before = nrow(input),
      n_locations_kept = sum(keep_location),
      n_locations_before = ncol(input),
      feature_counts = feature_counts,
      location_counts = location_counts,
      action = action
    )
    if (binary_skipped > 0L) {
      message(
        binary_skipped, " kept feature(s) have fewer than ", min_support,
        " detected or undetected locations; tweedspot(two_part = TRUE) will ",
        "test them with the count part only."
      )
    }
  }

  if (identical(action, "flag")) {
    output <- input
    rd <- SummarizedExperiment::rowData(output)
    rd$filter_keep <- keep_feature
    rd$filter_reason <- tweedspot_fail_reason(feature_fails)
    SummarizedExperiment::rowData(output) <- rd
    cd <- SummarizedExperiment::colData(output)
    cd$filter_keep <- keep_location
    cd$filter_reason <- tweedspot_fail_reason(location_fails)
    SummarizedExperiment::colData(output) <- cd
  } else {
    output <- input[keep_feature, keep_location]
  }
  if (!any(keep_feature) || !any(keep_location)) {
    warning("filter_features() kept ", sum(keep_feature), " feature(s) and ",
            sum(keep_location), " location(s); check the thresholds.", call. = FALSE)
  }
  S4Vectors::metadata(output)$filter_features <- summary_info
  output
}

tweedspot_location_mask <- function(input, location_mask) {
  if (is.null(location_mask)) {
    return(NULL)
  }
  cd <- SummarizedExperiment::colData(input)
  if (!is.character(location_mask) || length(location_mask) != 1L ||
      !location_mask %in% colnames(cd)) {
    stop("`location_mask` must name a column of `colData(input)`.")
  }
  mask <- cd[[location_mask]]
  if (is.numeric(mask) && all(mask[!is.na(mask)] %in% c(0, 1))) {
    mask <- mask == 1
  }
  if (!is.logical(mask)) {
    stop("`colData(input)$", location_mask, "` must be logical (or 0/1).")
  }
  mask[is.na(mask)] <- FALSE
  mask
}

# Per-group model support: a feature is supported in a group when it is
# detected in at least as many of that group's locations as the group's spatial
# smooth has basis functions. Groups below `min_group_size` are ignored.
tweedspot_group_support <- function(Y, coords, labels, smooth_k, min_group_size,
                                    detection_threshold = 0) {
  levels_ok <- levels(labels)[tabulate(labels, nlevels(labels)) >= min_group_size]
  supported <- matrix(FALSE, nrow(Y), length(levels_ok),
                      dimnames = list(rownames(Y), levels_ok))
  min_support <- stats::setNames(integer(length(levels_ok)), levels_ok)
  for (g in levels_ok) {
    idx <- which(labels == g)
    min_support[[g]] <- tweedspot_basis_k(coords[idx, , drop = FALSE], smooth_k)
    supported[, g] <- rowSums(Y[, idx, drop = FALSE] > detection_threshold) >= min_support[[g]]
  }
  list(supported = supported, min_support = min_support)
}
