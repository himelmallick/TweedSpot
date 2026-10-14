tweedspot_size_factor_methods <- function() {
  c("scran", "library_size", "spanorm")
}

tweedspot_libsize <- function(input, Y, assay_name = "counts",
                              size_factor_method = tweedspot_size_factor_methods(),
                              spanorm_args = list(),
                              spanorm_fit = NULL) {
  size_factor_method <- match.arg(size_factor_method)
  if (identical(size_factor_method, "spanorm")) {
    if (is.null(spanorm_fit)) {
      spanorm_fit <- tweedspot_run_spanorm(input, assay_name, spanorm_args)$fit
    }
    return(tweedspot_spanorm_offsets(
      spanorm_fit,
      gene_names = rownames(input),
      spot_names = colnames(input)
    ))
  }

  coldata <- SummarizedExperiment::colData(input)
  if ("sizeFactor" %in% names(coldata)) {
    sf <- tweedspot_as_numeric(coldata[["sizeFactor"]])
    if (length(sf) == ncol(input) && all(is.finite(sf) | is.na(sf))) {
      return(tweedspot_repair_size_factors(sf))
    }
  }

  if (identical(size_factor_method, "scran")) {
    sf <- tryCatch(
      scran::computeSumFactors(input, assay.type = assay_name),
      error = function(e) {
        stop(
          "SCRAN size factor estimation failed. This often happens when spatial locations ",
          "have zero library sizes; run `filter_features(..., drop_empty_locations = TRUE)`, ",
          "remove zero-count spots, or set `size_factor_method = \"library_size\"`. ",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    return(tweedspot_repair_size_factors(tweedspot_as_numeric(sf)))
  }
  tweedspot_library_size(Y)
}

# SpaNorm (Salim et al., Genome Biology 2025) fits, for every gene g and spot i,
#   Y_gi ~ NB(mu_gi, psi_g),
#   log mu_gi = gmean_g + f_g(x_i, y_i)                       [biology]
#             + {a_g + h_g(x_i, y_i)} log LS_i                 [library size]
#             + batch terms (optional)                         [batch]
# with thin-plate spline bases for f_g and h_g. It separates a spatially smooth
# library-size effect from spatial biology without domain or cluster labels.
#
# TweedSpot can use this fit in two ways (`spanorm_output`):
#   "size_factor": raw counts are modelled, and the technical part of the SpaNorm
#     linear predictor (library-size and batch columns) enters the GAM as a gene-
#     and spot-specific offset. The spatial biology is left for `s(x, y)` to test.
#   "counts": the SpaNorm percentile-adjusted counts (PAC) are modelled with a
#     unit offset. PAC are non-negative integers obtained by quantile-matching each
#     observed count from the full NB model to the biology-only NB model, so they
#     remain valid responses for every TweedSpot family.
tweedspot_spanorm_outputs <- function() {
  c("size_factor", "counts")
}

# A cached SpaNorm fit is reused only when it was made by TweedSpot for the same
# features, locations, assay values, SpaNorm settings, and seed.
tweedspot_spanorm_fingerprint <- function(input, assay_name, spanorm_args, seed) {
  Y <- SummarizedExperiment::assay(input, assay_name)
  list(
    features = rownames(input),
    locations = colnames(input),
    assay_name = assay_name,
    row_sums = as.numeric(rowSums(Y)),
    col_sums = as.numeric(colSums(Y)),
    weighted = as.numeric(colSums(Y * seq_len(nrow(Y)))),
    spanorm_args = spanorm_args,
    seed = seed
  )
}

tweedspot_run_spanorm <- function(input, assay_name = "counts",
                                  spanorm_args = list(), need_counts = FALSE,
                                  seed = NULL) {
  if (!requireNamespace("SpaNorm", quietly = TRUE)) {
    stop(
      "`size_factor_method = \"spanorm\"` requires the Bioconductor package `SpaNorm`. ",
      "Install it with `BiocManager::install(\"SpaNorm\")`.",
      call. = FALSE
    )
  }
  if (!is.list(spanorm_args)) {
    stop("`spanorm_args` must be a list of arguments passed to `SpaNorm::SpaNorm()`.")
  }
  if (!is.null(spanorm_args$adj.method) && !identical(spanorm_args$adj.method, "logpac")) {
    warning("`spanorm_args$adj.method` is ignored; TweedSpot uses SpaNorm's \"logpac\" adjustment.")
  }
  spanorm_args$adj.method <- NULL
  spanorm_args$spe <- NULL

  fingerprint <- tweedspot_spanorm_fingerprint(input, assay_name, spanorm_args, seed)
  restore_rng <- tweedspot_set_seed(seed)
  on.exit(restore_rng(), add = TRUE)
  md <- S4Vectors::metadata(input)
  reusable <- !is.null(md$SpaNorm) && identical(md$SpaNorm_fingerprint, fingerprint) &&
    tweedspot_spanorm_fit_matches(md$SpaNorm, input)
  if (reusable && !need_counts) {
    return(list(fit = md$SpaNorm, counts = NULL, fingerprint = fingerprint))
  }

  spe <- input
  if (!reusable) {
    # Never let SpaNorm silently reuse a fit made for other data or settings.
    S4Vectors::metadata(spe)$SpaNorm <- NULL
  }
  if (!identical(assay_name, "counts")) {
    SummarizedExperiment::assay(spe, "counts", withDimnames = FALSE) <-
      SummarizedExperiment::assay(input, assay_name, withDimnames = FALSE)
  }
  args <- c(list(spe = spe, adj.method = "logpac", verbose = FALSE), spanorm_args)
  spe <- tryCatch(
    suppressWarnings(do.call(SpaNorm::SpaNorm, args)),
    error = function(e) {
      stop("SpaNorm model fitting failed: ", conditionMessage(e),
           if (grepl("singular", conditionMessage(e))) {
             paste0(" This often happens with few locations; try ",
                    "`spanorm_args = list(sample.p = 1)` or a different `seed`.")
           },
           call. = FALSE)
    }
  )
  fit <- S4Vectors::metadata(spe)$SpaNorm
  if (!tweedspot_spanorm_fit_matches(fit, input)) {
    stop("Could not obtain a SpaNorm fit matching the genes and spots of `input`.")
  }
  counts <- NULL
  if (need_counts) {
    counts <- tweedspot_spanorm_pac(SummarizedExperiment::assay(spe, "logcounts"))
    dimnames(counts) <- dimnames(input)
  }
  list(fit = fit, counts = counts, fingerprint = fingerprint)
}

# SpaNorm's "logpac" output is log2(PAC + 1); PAC are integer NB quantiles.
tweedspot_spanorm_pac <- function(logpac) {
  pac <- round(2^as.matrix(logpac) - 1)
  pac[pac < 0] <- 0
  pac
}

tweedspot_spanorm_fit_matches <- function(fit, input) {
  !is.null(fit) &&
    identical(as.integer(fit$ngenes), nrow(input)) &&
    identical(as.integer(fit$ncells), ncol(input))
}

tweedspot_spanorm_offsets <- function(fit, gene_names = NULL, spot_names = NULL) {
  tech_cols <- which(fit$wtype %in% c("ls", "batch"))
  if (!any(fit$wtype == "ls")) {
    stop("The SpaNorm fit does not contain library-size terms.")
  }
  W_tech <- as.matrix(fit$W)[, tech_cols, drop = FALSE]
  alpha_tech <- as.matrix(fit$alpha)[, tech_cols, drop = FALSE]
  log_sf <- tcrossprod(alpha_tech, W_tech)
  # Centre each gene's log size factors; the constant is absorbed by the intercept.
  log_sf <- log_sf - rowMeans(log_sf)
  sf <- exp(log_sf)
  dimnames(sf) <- list(gene_names, spot_names)
  sf
}

tweedspot_as_numeric <- function(x) {
  if (methods::is(x, "SummarizedExperiment")) {
    coldata <- SummarizedExperiment::colData(x)
    if ("sizeFactor" %in% names(coldata)) {
      x <- coldata[["sizeFactor"]]
    } else {
      stop("Could not extract size factors from the object returned by `scran::computeSumFactors()`.")
    }
  }
  if (methods::is(x, "DataFrame") && ncol(x) == 1L) {
    x <- x[[1L]]
  }
  if (methods::is(x, "List") && length(x) == 1L) {
    x <- x[[1L]]
  }
  as.numeric(base::unlist(x, recursive = TRUE, use.names = FALSE))
}

tweedspot_repair_size_factors <- function(sf) {
  bad <- is.na(sf) | !is.finite(sf) | sf <= 0
  if (any(bad)) {
    warning(sum(bad), " size factor(s) were missing, zero, or negative and were set to 1. ",
            "Consider removing those locations (for example with filter_features()).",
            call. = FALSE)
  }
  sf[bad] <- 1
  sf
}

tweedspot_library_size <- function(Y) {
  ls <- colSums(Y)
  med_ls <- stats::median(ls)
  if (!is.finite(med_ls) || med_ls <= 0) {
    stop("Library sizes could not be computed because the median spot total is non-positive.")
  }
  ls / med_ls
}

tweedspot_check_scalar <- function(x, name, min_value = NULL, max_value = NULL) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single numeric value.", name))
  }
  if (!is.null(min_value) && x < min_value) {
    stop(sprintf("`%s` must be greater than or equal to %s.", name, min_value))
  }
  if (!is.null(max_value) && x > max_value) {
    stop(sprintf("`%s` must be less than or equal to %s.", name, max_value))
  }
}

tweedspot_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be `TRUE` or `FALSE`.", name))
  }
}

tweedspot_order_results <- function(res, key_col, decreasing = FALSE) {
  key <- res[[key_col]]
  stat <- res[["tweedspot_stat"]]
  if (decreasing) {
    ord <- order(is.na(key), -key, is.na(stat), -stat)
  } else {
    ord <- order(is.na(key), key, is.na(stat), -stat)
  }
  res[ord, , drop = FALSE]
}

# Readable-name columns searched (in order) when resolving and labelling features.
tweedspot_name_cols <- function(input) {
  intersect(c("feature_name", "gene_name"), colnames(SummarizedExperiment::rowData(input)))
}

tweedspot_resolve_gene <- function(input, gene) {
  rd <- SummarizedExperiment::rowData(input)
  name_cols <- tweedspot_name_cols(input)
  if (is.numeric(gene) && length(gene) == 1L && !is.na(gene)) {
    idx <- as.integer(gene)
    if (idx < 1L || idx > nrow(input)) {
      stop("Numeric `feature` index is out of bounds.")
    }
  } else if (is.character(gene) && length(gene) == 1L && !is.na(gene)) {
    idx <- match(gene, rownames(input))
    for (col in name_cols) {
      if (is.na(idx)) {
        idx <- match(gene, as.character(rd[[col]]))
      }
    }
    if (is.na(idx)) {
      stop("`feature` was not found among row names or `rowData(input)` ",
           "feature_name/gene_name columns.")
    }
  } else {
    stop("`feature` must be a single row index, row name, or feature name.")
  }

  id <- rownames(input)[idx]
  if (is.null(id)) {
    id <- paste0("feature", idx)
  }
  label <- id
  for (col in rev(name_cols)) {
    nm <- as.character(rd[[col]][idx])
    if (!is.na(nm) && nzchar(nm)) {
      label <- if (identical(nm, id)) id else paste0(nm, " (", id, ")")
    }
  }
  list(index = idx, label = label)
}

tweedspot_feature_names <- function(input, feature_name_col = NULL) {
  if (is.null(feature_name_col)) {
    nm <- rownames(input)
    if (is.null(nm)) {
      nm <- rep("", nrow(input))
    }
  } else {
    rd <- SummarizedExperiment::rowData(input)
    if (!feature_name_col %in% colnames(rd)) {
      stop(sprintf("`feature_name_col` not found in `rowData(input)`: %s.", feature_name_col))
    }
    nm <- as.character(rd[[feature_name_col]])
  }
  nm[is.na(nm)] <- ""
  nm
}

# Modality-agnostic exclusions: an explicit list of feature names and/or regular
# expressions, matched against row names and, if given, `feature_name_col`.
tweedspot_explicit_exclusions <- function(input, exclude = NULL, exclude_pattern = NULL,
                                          feature_name_col = NULL) {
  drop <- rep(FALSE, nrow(input))
  if (is.null(exclude) && is.null(exclude_pattern)) {
    return(drop)
  }
  ids <- rownames(input)
  if (is.null(ids)) {
    ids <- rep("", nrow(input))
  }
  labels <- list(ids)
  if (!is.null(feature_name_col)) {
    labels <- c(labels, list(tweedspot_feature_names(input, feature_name_col)))
  }
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      stop("`exclude` must be a character vector of feature names.")
    }
    for (lab in labels) {
      drop <- drop | lab %in% exclude
    }
  }
  if (!is.null(exclude_pattern)) {
    if (!is.character(exclude_pattern) || anyNA(exclude_pattern)) {
      stop("`exclude_pattern` must be a character vector of regular expressions.")
    }
    for (pat in exclude_pattern) {
      for (lab in labels) {
        drop <- drop | grepl(pat, lab)
      }
    }
  }
  drop
}

# Transcriptomics-specific exclusions for mitochondrial and ribosomal genes.
tweedspot_annotation_exclusions <- function(input, exclude_mito, exclude_ribo,
                                            feature_name_col = NULL,
                                            feature_type_col = NULL) {
  feature_names <- tweedspot_feature_names(input, feature_name_col)
  types <- NULL
  if (!is.null(feature_type_col)) {
    rd <- SummarizedExperiment::rowData(input)
    if (!feature_type_col %in% colnames(rd)) {
      stop(sprintf("`feature_type_col` not found in `rowData(input)`: %s.", feature_type_col))
    }
    types <- as.character(rd[[feature_type_col]])
    types[is.na(types)] <- ""
  }

  none <- rep(FALSE, nrow(input))
  mito <- none
  ribo <- none

  if (exclude_mito) {
    # Human "MT-", mouse "mt-", and "MT." after make.names(); not e.g. MTOR.
    mito_name <- grepl("^mt[-.]", feature_names, ignore.case = TRUE)
    mito_type <- if (is.null(types)) none else grepl("mitochond", types, ignore.case = TRUE)
    mito <- mito_name | mito_type
  }

  if (exclude_ribo) {
    # RPL/RPS/MRPL/MRPS protein genes (RPL10A, RPS4Y1, RPLP0, RPSA, MRPL12, Rpl7a)
    # but not look-alikes such as the RPS6KA/RPS6KB kinases.
    ribo_name <- grepl(
      "^(rpl|rps|mrpl|mrps)([0-9]+[a-z]?[0-9]*|p[0-9]|a)$",
      feature_names,
      ignore.case = TRUE
    )
    ribo_type <- if (is.null(types)) none else grepl("ribosom", types, ignore.case = TRUE)
    ribo <- ribo_name | ribo_type
  }

  mito | ribo
}

tweedspot_feature_filters <- function(Y, min_total, min_pct_locations,
                                      min_locations, min_mean, min_var,
                                      detection_threshold = 0,
                                      min_support = 0L) {
  detected <- Y > detection_threshold
  n_detected <- rowSums(detected)
  fails <- list(
    low_total = rowSums(Y) < min_total,
    low_pct_locations = rowMeans(detected) * 100 < min_pct_locations,
    low_locations = n_detected < min_locations,
    low_mean = rowMeans(Y) < min_mean,
    low_var = tweedspot_row_vars(Y) < min_var,
    low_model_support = n_detected < min_support
  )
  lapply(fails, function(x) {
    x[is.na(x)] <- TRUE
    x
  })
}

tweedspot_location_filters <- function(Y, detection_threshold = 0,
                                       min_total_per_location = 0,
                                       min_features_per_location = 0,
                                       outlier_nmads = NULL,
                                       excluded_fraction = NULL,
                                       max_excluded_fraction = NULL,
                                       mask = NULL) {
  total <- colSums(Y)
  n_features <- colSums(Y > detection_threshold)
  none <- rep(FALSE, ncol(Y))
  fails <- list(
    low_total = total < min_total_per_location,
    low_features = n_features < min_features_per_location,
    outlier = none,
    high_excluded_fraction = none,
    masked = none
  )
  if (!is.null(outlier_nmads)) {
    fails$outlier <- tweedspot_low_outlier(log1p(total), outlier_nmads) |
      tweedspot_low_outlier(log1p(n_features), outlier_nmads)
  }
  if (!is.null(max_excluded_fraction) && !is.null(excluded_fraction)) {
    fails$high_excluded_fraction <- excluded_fraction > max_excluded_fraction
  }
  if (!is.null(mask)) {
    fails$masked <- !mask
  }
  lapply(fails, function(x) {
    x[is.na(x)] <- TRUE
    x
  })
}

# Lower-tail outliers: more than `nmads` median absolute deviations below the
# median. Only the lower tail is used, as low signal indicates poor quality in
# every modality while high signal can be genuine biology.
tweedspot_low_outlier <- function(x, nmads) {
  med <- stats::median(x)
  spread <- stats::mad(x)
  if (!is.finite(spread) || spread == 0) {
    return(rep(FALSE, length(x)))
  }
  x < med - nmads * spread
}

tweedspot_fail_reason <- function(fails) {
  M <- do.call(cbind, fails)
  if (is.null(dim(M))) {
    M <- matrix(M, nrow = 1L, dimnames = list(NULL, names(fails)))
  }
  unname(apply(M, 1, function(row) paste(colnames(M)[row], collapse = ";")))
}

tweedspot_report_filtering <- function(n_features_kept, n_features_before,
                                       n_locations_kept, n_locations_before,
                                       feature_counts, location_counts,
                                       action) {
  verb <- if (identical(action, "flag")) "flagged to keep" else "retained"
  message(
    "filter_features ", verb, " ", n_features_kept, "/", n_features_before,
    " features and ", n_locations_kept, "/", n_locations_before, " spatial locations."
  )
  fmt <- function(x) {
    x <- x[x > 0]
    if (!length(x)) "none" else paste(names(x), x, sep = "=", collapse = ", ")
  }
  message("Features failing each criterion: ", fmt(feature_counts), ".")
  message("Locations failing each criterion: ", fmt(location_counts), ".")
}

# Parse `covariates` into the variables it uses and the right-hand-side
# expression. Character vectors are combined with `+`; names are kept exactly,
# so columns such as `conc (mM)` can be used.
tweedspot_parse_covariates <- function(covariates) {
  if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L) {
      stop("`covariates` must be a one-sided formula like `~ batch + age`.")
    }
    rhs <- covariates[[2L]]
  } else if (is.character(covariates) && length(covariates) && !anyNA(covariates)) {
    rhs <- Reduce(function(a, b) call("+", a, b), lapply(covariates, as.name))
  } else {
    stop("`covariates` must be NULL, a one-sided formula, or a character vector of `colData(input)` column names.")
  }
  list(vars = all.vars(rhs), rhs = rhs)
}

# Replace variable names inside an expression, leaving function names alone.
tweedspot_rename_symbols <- function(expr, map) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    return(if (nm %in% names(map)) as.name(map[[nm]]) else expr)
  }
  if (is.call(expr) && length(expr) > 1L) {
    for (i in seq_along(expr)[-1L]) {
      if (!is.null(expr[[i]])) {
        expr[[i]] <- tweedspot_rename_symbols(expr[[i]], map)
      }
    }
  }
  expr
}

tweedspot_validate_covariate_data <- function(covariate_data) {
  bad_type <- vapply(covariate_data, function(x) {
    !(is.numeric(x) || is.logical(x) || is.character(x) || is.factor(x))
  }, logical(1))
  if (any(bad_type)) {
    stop("Referenced `covariates` columns must be numeric, logical, character, or factor.")
  }
  if (anyNA(covariate_data)) {
    stop("Referenced `covariates` columns must not contain missing values.")
  }
  for (j in seq_along(covariate_data)) {
    if (is.character(covariate_data[[j]])) {
      covariate_data[[j]] <- as.factor(covariate_data[[j]])
    }
  }
  constant_covariates <- vapply(covariate_data, function(x) {
    length(unique(x)) < 2L
  }, logical(1))
  if (any(constant_covariates)) {
    stop(sprintf(
      "Referenced `covariates` must vary across spatial locations. Constant columns: %s.",
      paste(names(covariate_data)[constant_covariates], collapse = ", ")
    ))
  }
  covariate_data
}

tweedspot_sanitize_covariates <- function(covariate_data, covariate_vars, rhs) {
  reserved <- c("expr", "x", "y", "libsz")
  sanitized_names <- make.names(covariate_vars, unique = TRUE)
  if (any(sanitized_names %in% reserved)) {
    stop("Referenced `covariates` columns must not use reserved names: expr, x, y, libsz.")
  }
  names(covariate_data) <- sanitized_names
  map <- stats::setNames(sanitized_names, covariate_vars)
  terms <- paste(deparse(tweedspot_rename_symbols(rhs, map), width.cutoff = 500L),
                 collapse = " ")
  list(data = covariate_data, terms = terms)
}

tweedspot_covariates <- function(input, covariates) {
  if (is.null(covariates)) {
    return(NULL)
  }

  # `optional = TRUE` keeps column names exactly as they are in colData().
  coldata <- as.data.frame(SummarizedExperiment::colData(input),
                           optional = TRUE, stringsAsFactors = FALSE)
  parsed <- tweedspot_parse_covariates(covariates)
  covariate_vars <- parsed$vars
  if (!length(covariate_vars)) {
    stop("`covariates` must reference at least one column in `colData(input)`.")
  }

  missing_vars <- setdiff(covariate_vars, colnames(coldata))
  if (length(missing_vars)) {
    stop(sprintf("`covariates` references columns not found in `colData(input)`: %s.",
                 paste(missing_vars, collapse = ", ")))
  }

  covariate_data <- coldata[, covariate_vars, drop = FALSE]
  covariate_data <- tweedspot_validate_covariate_data(covariate_data)
  tweedspot_sanitize_covariates(covariate_data, covariate_vars, parsed$rhs)
}

tweedspot_resolve_family <- function(family, tweedie_p = NULL, Y = NULL) {
  if (is.null(tweedie_p) || length(tweedie_p) == 0L ||
      (length(tweedie_p) == 1L && is.na(tweedie_p))) {
    if (identical(family, "tw")) {
      # Same bound as `tweedie_p = "auto"`: on integer data the Tweedie density
      # spikes at the integers as p -> 1, which inflates the smooth's EDF and
      # makes the spatial test anti-conservative.
      a <- if (!is.null(Y) && tweedspot_is_integer_valued(Y)) {
        tweedspot_auto_min_p()
      } else {
        1.01
      }
      family <- tweedspot_tw(a = a)
    }
    return(tweedspot_mark_family(family))
  }
  if (is.character(tweedie_p)) {
    if (length(tweedie_p) == 1L && identical(tolower(tweedie_p), "auto")) {
      return("auto")
    }
    stop("`tweedie_p` must be NULL, \"auto\", or one of 0, 1, values in (1, 2], or 3.")
  }
  tweedspot_check_scalar(tweedie_p, "tweedie_p")
  if (tweedie_p > 0 && tweedie_p < 1) {
    stop("`tweedie_p` values between 0 and 1 are not supported.")
  }
  if (tweedie_p < 0 || tweedie_p > 3 || (tweedie_p > 2 && tweedie_p < 3)) {
    stop("`tweedie_p` must be one of 0, 1, values in (1, 2], or 3.")
  }
  if (tweedie_p == 0) {
    # Log link so that offset(log(libsz)) acts multiplicatively, as for every
    # other family supported by TweedSpot.
    return(stats::gaussian(link = "log"))
  }
  if (tweedie_p == 1) {
    return(stats::poisson(link = "log"))
  }
  if (tweedie_p > 1 && tweedie_p <= 2) {
    if (tweedie_p == 2 && !is.null(Y) && min(Y, na.rm = TRUE) <= 0) {
      stop("`tweedie_p = 2` corresponds to a Gamma model and requires strictly positive responses.")
    }
    return(mgcv::Tweedie(p = tweedie_p, link = "log"))
  }
  if (tweedie_p == 3) {
    if (!is.null(Y) && min(Y, na.rm = TRUE) <= 0) {
      stop("`tweedie_p = 3` uses inverse Gaussian models and requires strictly positive responses.")
    }
    return(stats::inverse.gaussian(link = "log"))
  }
}

# `mgcv::tw()` evaluates helpers such as `ldTweedie()` from an environment whose
# parent is the global environment, so fits fail unless mgcv is attached (and
# always fail on BiocParallel workers). Re-parent it to the mgcv namespace.
tweedspot_tw <- function(link = "log", a = 1.01, b = 1.99) {
  fam <- mgcv::tw(link = link, a = a, b = b)
  env <- environment(fam$ls)
  if (!is.null(env) && identical(parent.env(env), globalenv())) {
    parent.env(env) <- asNamespace("mgcv")
  }
  fam
}

# Candidate variance powers for `tweedie_p = "auto"`. All candidates share the
# log link, the spatial smooth, covariates, and offset, so they differ only in
# the variance function V(mu) = phi * mu^p, and every candidate estimates its
# own dispersion phi. Powers 2 (Gamma) and 3 (inverse Gaussian) have no mass at
# zero and are only considered for strictly positive responses.
#
# Poisson (p = 1, phi fixed at 1) is deliberately not a candidate: because the
# family is chosen with the spatial smooth in the model, genes whose noise
# happens to look spatial fit Poisson unusually well, get selected as Poisson,
# and are then tested with no allowance for overdispersion, which inflated
# false positives (about 10% at alpha = 0.05 among those genes). Poisson
# remains available as a fixed choice, `tweedie_p = 1`.
#
# For integer data the estimated-Tweedie candidate is restricted to p >= 1.1:
# as p -> 1 with phi near 1, the Tweedie density on integer data collapses onto
# spikes at the integers, which inflates its likelihood and makes the spatial
# test anti-conservative.
tweedspot_auto_candidates <- function(expr) {
  cands <- list(gaussian = stats::gaussian(link = "log"))
  cands$tweedie <- if (tweedspot_is_integer_valued(expr)) {
    tweedspot_tw(a = tweedspot_auto_min_p())
  } else {
    tweedspot_tw()
  }
  if (all(expr > 0, na.rm = TRUE)) {
    cands$gamma <- stats::Gamma(link = "log")
    cands$inverse_gaussian <- stats::inverse.gaussian(link = "log")
  }
  cands
}

tweedspot_auto_min_p <- function() 1.1

# Works for base vectors and matrices, Matrix sparse matrices, and DelayedArray
# objects without making a dense copy of the whole assay.
tweedspot_is_integer_valued <- function(x, tol = 1e-8) {
  if (methods::is(x, "sparseMatrix")) {
    x <- methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix")@x
  } else if (methods::is(x, "DelayedArray")) {
    if (!length(x)) {
      return(TRUE)
    }
    return(max(abs(x - round(x)), na.rm = TRUE) < tol)
  }
  x <- x[!is.na(x)]
  all(abs(x - round(x)) < tol)
}

# The assay to model, kept in its own storage (dense, sparse, or DelayedArray)
# so large data are never copied into one dense matrix.
tweedspot_get_assay <- function(input, assay_name) {
  if (!assay_name %in% SummarizedExperiment::assayNames(input)) {
    stop("Assay `", assay_name, "` not found in `input`.", call. = FALSE)
  }
  Y <- SummarizedExperiment::assay(input, assay_name)
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }
  if (anyNA(Y)) {
    stop("Assay `", assay_name, "` contains missing values; TweedSpot needs a complete ",
         "feature-by-location matrix.", call. = FALSE)
  }
  Y
}

tweedspot_row_vars <- function(Y) {
  n <- ncol(Y)
  if (n < 2L) {
    return(rep(NA_real_, nrow(Y)))
  }
  s1 <- rowSums(Y)
  unname(pmax((rowSums(Y^2) - s1^2 / n) / (n - 1), 0))
}

tweedspot_family_power <- function(fit) {
  fam <- fit$family
  name <- fam$family
  if (is.null(name)) {
    return(NA_real_)
  }
  if (identical(name, "gaussian")) return(0)
  if (name %in% c("poisson", "quasipoisson")) return(1)
  if (identical(name, "Gamma")) return(2)
  if (identical(name, "inverse.gaussian")) return(3)
  if (grepl("^Tweedie", name)) {
    if (is.function(fam$getTheta)) {
      return(as.numeric(fam$getTheta(TRUE)))
    }
    p <- tryCatch(get("p", envir = environment(fam$variance)), error = function(e) NA_real_)
    return(as.numeric(p))
  }
  NA_real_
}

tweedspot_information_criterion <- function(fit, criterion = c("BIC", "AIC")) {
  criterion <- match.arg(criterion)
  ll <- stats::logLik(fit)
  df <- attr(ll, "df")
  n <- length(fit$y)
  value <- as.numeric(ll)
  if (!is.finite(value)) {
    return(Inf)
  }
  penalty <- if (identical(criterion, "BIC")) log(n) else 2
  -2 * value + penalty * df
}

# Extended families such as tw() and nb() keep their current parameter (the
# Tweedie power or NB theta) in a closure, and mgcv starts each fit from it. A
# family object reused across genes therefore carries one gene's estimate into
# the next fit, so results depend on gene order and can get stuck at a boundary.
# Record the starting value once and restore it before every fit.
tweedspot_mark_family <- function(family) {
  if (inherits(family, "extended.family") && is.function(family$getTheta) &&
      is.function(family$putTheta)) {
    attr(family, "tweedspot_theta0") <- family$getTheta()
  }
  family
}

tweedspot_reset_family <- function(family) {
  theta0 <- attr(family, "tweedspot_theta0")
  if (!is.null(theta0)) {
    family$putTheta(theta0)
  }
  family
}

tweedspot_fit <- function(formula, family, data, fit_method, use_bam,
                          bam_discrete, bam_nthreads, p_criterion = "BIC") {
  family <- tweedspot_reset_family(family)
  if (identical(family, "auto")) {
    return(tweedspot_fit_auto(formula, data, fit_method, use_bam,
                              bam_discrete, bam_nthreads, p_criterion))
  }
  fit_once <- function(...) {
    if (use_bam) {
      return(mgcv::bam(
        formula = formula,
        family = family,
        data = data,
        method = fit_method,
        discrete = bam_discrete,
        nthreads = bam_nthreads,
        ...
      ))
    }
    mgcv::gam(
      formula = formula,
      family = family,
      data = data,
      method = fit_method,
      ...
    )
  }
  if (tweedspot_is_log_gaussian(family)) {
    # Gaussian with log link can fail to find starting values on sparse counts.
    return(tryCatch(fit_once(), error = function(e) {
      y <- data[[all.vars(formula)[1L]]]
      fit_once(mustart = pmax(y, mean(y) / 10, 1e-3))
    }))
  }
  fit_once()
}

tweedspot_is_log_gaussian <- function(family) {
  inherits(family, "family") && identical(family$family, "gaussian") &&
    identical(family$link, "log")
}

# Fit each candidate variance power and keep the one with the smallest
# information criterion. Returns the selected fit with attributes recording the
# criterion value of every candidate.
tweedspot_fit_auto <- function(formula, data, fit_method, use_bam, bam_discrete,
                               bam_nthreads, p_criterion = "BIC") {
  expr <- data[[all.vars(formula)[1L]]]
  cands <- tweedspot_auto_candidates(expr)
  fits <- lapply(cands, function(fam) {
    tryCatch(
      suppressWarnings(tweedspot_fit(formula, fam, data, fit_method, use_bam,
                                     bam_discrete, bam_nthreads)),
      error = function(e) NULL
    )
  })
  ic <- vapply(fits, function(f) {
    if (is.null(f)) Inf else tweedspot_information_criterion(f, p_criterion)
  }, numeric(1))
  if (!any(is.finite(ic))) {
    stop("All candidate families failed for `tweedie_p = \"auto\"`.")
  }
  best <- which.min(ic)
  fit <- fits[[best]]
  attr(fit, "tweedspot_ic") <- ic
  attr(fit, "tweedspot_selected") <- names(cands)[best]
  fit
}

tweedspot_fit_single_gene <- function(expr, coords, libsz, covariates, family,
                                      fit_method, use_bam, bam_discrete,
                                      bam_nthreads, smooth_k, p_criterion = "BIC") {
  tweedspot_validate_libsize(libsz)
  formula_parts <- tweedspot_model_formulas(coords, libsz, covariates, smooth_k)
  form <- formula_parts$form
  d_base <- formula_parts$data
  d <- cbind(expr = as.numeric(expr), d_base)
  tryCatch(
    tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete, bam_nthreads,
                  p_criterion = p_criterion),
    error = function(e) NULL
  )
}

tweedspot_basis_k <- function(coords, smooth_k = NULL, default_k = 30L) {
  n_unique <- nrow(unique(as.data.frame(coords[, seq_len(min(2, ncol(coords))), drop = FALSE])))
  if (n_unique <= 3L) {
    return(3L)
  }
  if (is.null(smooth_k)) {
    return(as.integer(min(default_k, n_unique - 1L)))
  }
  as.integer(max(3L, min(as.integer(smooth_k), n_unique - 1L)))
}

tweedspot_validate_libsize <- function(libsz) {
  if (any(!is.finite(libsz) | libsz <= 0)) {
    stop(
      "Library sizes must be positive and finite for every spatial location. ",
      "This often happens when filtered spots have zero total counts; remove those spots ",
      "or provide a positive `colData(input)$sizeFactor`."
    )
  }
}

tweedspot_model_frame <- function(coords, libsz, covariates) {
  # Gene-specific size factors (a genes x spots matrix) are inserted per gene in
  # `tweedspot_fit_row()`; the frame carries a placeholder until then.
  libsz_col <- if (is.matrix(libsz)) rep(1, nrow(coords)) else libsz
  d_base <- data.frame(x = coords[, 1], y = coords[, 2], libsz = libsz_col)
  if (!is.null(covariates)) {
    d_base <- cbind(d_base, covariates$data)
  }
  d_base
}

tweedspot_model_formulas <- function(coords, libsz, covariates, smooth_k) {
  d_base <- tweedspot_model_frame(coords, libsz, covariates)
  basis_k <- tweedspot_basis_k(coords, smooth_k)
  smooth_term <- sprintf("s(x, y, k = %d)", basis_k)
  cov_terms <- if (is.null(covariates)) character(0) else covariates$terms
  rhs <- paste(c(smooth_term, cov_terms, "offset(log(libsz))"), collapse = " + ")
  list(
    data = d_base,
    form = stats::as.formula(sprintf("expr ~ %s", rhs)),
    form_bin = stats::as.formula(sprintf("I(expr > 0) ~ %s", rhs))
  )
}

# The presence/absence part of the two-part model needs at least `min_support`
# detected and `min_support` undetected locations to estimate a spatial smooth
# with that many basis functions; otherwise only the count part is used.
tweedspot_binary_estimable <- function(expr, min_support) {
  expr <- expr[!is.na(expr)]
  sum(expr > 0) >= min_support && sum(expr <= 0) >= min_support
}

tweedspot_smooth_summary <- function(fit, two_part, form_bin, data, fit_method,
                                     use_bam, bam_discrete, bam_nthreads,
                                     min_support = 0L) {
  s <- summary(fit)$s.table
  model_summary <- summary(fit)
  plog <- NA_real_
  if (two_part && tweedspot_binary_estimable(data$expr, min_support)) {
    f1 <- tryCatch(
      tweedspot_fit(form_bin, "binomial", data, fit_method, use_bam,
                    bam_discrete, bam_nthreads),
      error = function(e) NULL
    )
    if (!is.null(f1)) {
      plog <- summary(f1)$s.table[1, "p-value"]
    }
  }
  c(
    stat = s[1, 3],
    pval = s[1, "p-value"],
    plogit = plog,
    edf = s[1, 1],
    dev_expl = model_summary$dev.expl,
    p = tweedspot_family_power(fit)
  )
}

tweedspot_empty_row <- function() {
  c(
    stat = NA_real_,
    pval = NA_real_,
    plogit = NA_real_,
    edf = NA_real_,
    dev_expl = NA_real_,
    p = NA_real_
  )
}

tweedspot_fit_row <- function(i, Y, d_base, form, form_bin, two_part, family,
                              fit_method, use_bam, bam_discrete, bam_nthreads,
                              libsz = NULL, p_criterion = "BIC", min_support = 0L) {
  d <- cbind(expr = as.numeric(Y[i, ]), d_base)
  if (is.matrix(libsz)) {
    d$libsz <- as.numeric(libsz[i, ])
  }
  fit <- tryCatch(
    tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete, bam_nthreads,
                  p_criterion = p_criterion),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tweedspot_empty_row())
  }
  tweedspot_smooth_summary(
    fit = fit,
    two_part = two_part,
    form_bin = form_bin,
    data = d,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    bam_nthreads = bam_nthreads,
    min_support = min_support
  )
}

tweedspot_collect_fits <- function(Y, d_base, form, form_bin, two_part, family,
                                   fit_method, use_bam, bam_discrete,
                                   bam_nthreads, BPPARAM, libsz = NULL,
                                   p_criterion = "BIC", min_support = 0L) {
  if (methods::is(Y, "CsparseMatrix")) {
    # Row-compressed storage makes extracting one feature at a time fast.
    Y <- methods::as(Y, "RsparseMatrix")
  }
  one <- function(i) {
    tweedspot_fit_row(
      i = i,
      Y = Y,
      d_base = d_base,
      form = form,
      form_bin = form_bin,
      two_part = two_part,
      family = family,
      fit_method = fit_method,
      use_bam = use_bam,
      bam_discrete = bam_discrete,
      bam_nthreads = bam_nthreads,
      libsz = libsz,
      p_criterion = p_criterion,
      min_support = min_support
    )
  }
  do.call(rbind, BiocParallel::bplapply(seq_len(nrow(Y)), one, BPPARAM = BPPARAM))
}

tweedspot_extract_one_part <- function(M) {
  if (all(is.na(M[, "pval"]))) {
    stop(
      "All gene-wise TweedSpot fits failed. Check covariates, library sizes, and filtering thresholds."
    )
  }
  list(
    stat = M[, "stat"],
    pval = M[, "pval"],
    edf = M[, "edf"],
    dev_expl = M[, "dev_expl"],
    p = M[, "p"]
  )
}

# Two-part p-value: Cauchy combination of the count and presence/absence parts.
tweedspot_combine_pvalues <- function(M) {
  vapply(seq_len(nrow(M)), function(i) {
    if (is.na(M[i, "pval"])) {
      return(NA_real_)
    }
    # A skipped or failed presence/absence part contributes no evidence.
    cct(stats::na.omit(c(M[i, "pval"], M[i, "plogit"])))
  }, numeric(1))
}

tweedspot_extract_two_part <- function(M) {
  combined <- tweedspot_combine_pvalues(M)
  if (all(is.na(combined))) {
    stop(
      "All gene-wise TweedSpot fits failed. Check covariates, library sizes, and filtering thresholds."
    )
  }
  list(
    stat = M[, "stat"],
    pval = combined,
    edf = M[, "edf"],
    dev_expl = M[, "dev_expl"],
    p = M[, "p"]
  )
}

tweedspot_agnostic <- function(Y, coords, libsz, covariates, two_part,
                               family, fit_method, use_bam, bam_discrete,
                               bam_nthreads, smooth_k, BPPARAM,
                               p_criterion = "BIC") {
  tweedspot_validate_libsize(libsz)
  if (is.matrix(libsz) && !identical(dim(libsz), dim(Y))) {
    stop("Gene-specific size factors must be a matrix with the same dimensions as the counts.")
  }
  formula_parts <- tweedspot_model_formulas(coords, libsz, covariates, smooth_k)
  M <- tweedspot_collect_fits(
    Y = Y,
    d_base = formula_parts$data,
    form = formula_parts$form,
    form_bin = formula_parts$form_bin,
    two_part = two_part,
    family = family,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    bam_nthreads = bam_nthreads,
    BPPARAM = BPPARAM,
    libsz = libsz,
    p_criterion = p_criterion,
    min_support = tweedspot_basis_k(coords, smooth_k)
  )
  if (!two_part) {
    return(tweedspot_extract_one_part(M))
  }
  tweedspot_extract_two_part(M)
}

# Spatial coordinates for model fitting: the first two columns, centred and
# divided by ONE common scale. A common scale keeps the tissue's geometry, which
# matters because the thin-plate smooth s(x, y) is isotropic; scaling each axis
# separately would stretch elongated tissues and cost power.
tweedspot_model_coords <- function(input) {
  coords <- SpatialExperiment::spatialCoords(input)
  if (ncol(coords) < 2L) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  if (ncol(coords) > 2L) {
    warning("`spatialCoords(input)` has ", ncol(coords), " columns; TweedSpot uses the ",
            "first two (", paste(colnames(coords)[1:2], collapse = ", "), ").", call. = FALSE)
  }
  xy <- coords[, 1:2, drop = FALSE]
  storage.mode(xy) <- "double"
  missing <- !stats::complete.cases(xy)
  if (any(missing)) {
    stop(sum(missing), " location(s) have missing spatial coordinates; remove them before ",
         "running TweedSpot.", call. = FALSE)
  }
  spread <- apply(xy, 2, stats::sd)
  if (nrow(xy) > 1L && any(spread == 0)) {
    stop("Spatial coordinate `", colnames(xy)[spread == 0][1L], "` is constant; ",
         "TweedSpot needs locations that vary in two dimensions.", call. = FALSE)
  }
  common <- sqrt(mean(spread^2))
  if (!is.finite(common) || common == 0) {
    common <- 1
  }
  xy <- sweep(xy, 2, colMeans(xy)) / common
  colnames(xy) <- c("x", "y")
  xy
}

tweedspot_prepare_input <- function(input, assay_name, covariates,
                                    size_factor_method, spanorm_args = list(),
                                    spanorm_output = "size_factor", verbose = FALSE,
                                    seed = NULL) {
  Y <- tweedspot_get_assay(input, assay_name)
  coords <- tweedspot_model_coords(input)
  spanorm_fit <- NULL
  spanorm_fingerprint <- NULL
  if (identical(size_factor_method, "spanorm")) {
    spanorm_output <- match.arg(spanorm_output, tweedspot_spanorm_outputs())
    if (verbose) {
      message(
        "Note: SpaNorm's library-size effect is estimated from the same data and ",
        "treated as known. If library size varies across the tissue, spatial tests ",
        "can be anti-conservative; check calibration (see ?tweedspot)."
      )
    }
    sn <- tweedspot_run_spanorm(
      input = input,
      assay_name = assay_name,
      spanorm_args = spanorm_args,
      need_counts = identical(spanorm_output, "counts"),
      seed = seed
    )
    spanorm_fit <- sn$fit
    spanorm_fingerprint <- sn$fingerprint
    if (identical(spanorm_output, "counts")) {
      # PAC already have library-size effects removed, so no further offset.
      Y <- sn$counts
      libsz <- rep(1, ncol(Y))
    } else {
      libsz <- tweedspot_spanorm_offsets(
        spanorm_fit,
        gene_names = rownames(input),
        spot_names = colnames(input)
      )
    }
  } else {
    libsz <- tweedspot_libsize(
      input = input,
      Y = Y,
      assay_name = assay_name,
      size_factor_method = size_factor_method
    )
  }
  list(
    Y = Y,
    coords = coords,
    libsz = libsz,
    covariates = tweedspot_covariates(input, covariates),
    spanorm_fit = spanorm_fit,
    spanorm_fingerprint = spanorm_fingerprint
  )
}

# Set the random seed for the duration of a call and return a function that
# restores the caller's random-number state afterwards.
tweedspot_set_seed <- function(seed) {
  if (is.null(seed)) {
    return(function() invisible(NULL))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != round(seed)) {
    stop("`seed` must be NULL or a single whole number.")
  }
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}

# Remove every result a previous tweedspot() run wrote, so each run starts clean.
tweedspot_clear_results <- function(input, verbose = FALSE) {
  rd <- SummarizedExperiment::rowData(input)
  old_cols <- grepl("^tweedspot_", colnames(rd))
  md <- S4Vectors::metadata(input)
  had_results <- any(old_cols) || !is.null(md$tweedspot) || !is.null(md$tweedspot_group)
  if (any(old_cols)) {
    SummarizedExperiment::rowData(input) <- rd[, !old_cols, drop = FALSE]
  }
  md$tweedspot <- NULL
  md$tweedspot_group <- NULL
  S4Vectors::metadata(input) <- md
  if (had_results && verbose) {
    message("Replacing results from a previous tweedspot() run.")
  }
  input
}

# Warn when locations were subset or reordered after the run that produced the
# stored results: those results still describe the original locations.
tweedspot_check_locations <- function(input, stored_locations) {
  if (is.null(stored_locations)) {
    return(invisible(TRUE))
  }
  if (!identical(colnames(input), stored_locations)) {
    warning("Spatial locations were subset or reordered after `tweedspot()`; ",
            "the stored results describe the original ", length(stored_locations),
            " locations. Re-run `tweedspot()` on the current object.", call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

tweedspot_message_start <- function(verbose, Y, BPPARAM) {
  if (verbose) {
    message(
      "Running TweedSpot on ", nrow(Y), " genes across ", ncol(Y), " spatial locations ",
      "with ", BiocParallel::bpnworkers(BPPARAM), " worker(s)."
    )
  }
}

tweedspot_store_results <- function(input, res, padj_method) {
  SummarizedExperiment::rowData(input)$tweedspot_stat <- res$stat
  SummarizedExperiment::rowData(input)$tweedspot_pval <- res$pval
  SummarizedExperiment::rowData(input)$tweedspot_padj <-
    stats::p.adjust(res$pval, method = padj_method)
  SummarizedExperiment::rowData(input)$tweedspot_edf <- res$edf
  SummarizedExperiment::rowData(input)$tweedspot_dev_expl <- res$dev_expl
  SummarizedExperiment::rowData(input)$tweedspot_p <- res$p
  SummarizedExperiment::rowData(input)$tweedspot_family <- tweedspot_power_label(res$p)
  input
}

tweedspot_power_label <- function(p) {
  out <- rep(NA_character_, length(p))
  ok <- !is.na(p)
  out[ok & p == 0] <- "gaussian"
  out[ok & p == 1] <- "poisson"
  out[ok & p > 1 & p < 2] <- "tweedie"
  out[ok & p == 2] <- "gamma"
  out[ok & p == 3] <- "inverse_gaussian"
  out
}

cct <- function(pvals, weights = NULL) {
  pvals <- ifelse(is.na(pvals), 1, pvals)
  pvals <- pmin(pmax(pvals, .Machine$double.xmin), 1 - .Machine$double.eps)

  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else {
    weights <- weights / sum(weights)
  }

  is.small <- pvals < 1e-16
  cct <- if (!any(is.small)) {
    sum(weights * tan((0.5 - pvals) * pi))
  } else {
    sum((weights[is.small] / pvals[is.small]) / pi) +
      sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  if (cct > 1e15) (1 / cct) / pi else 1 - stats::pcauchy(cct)
}
