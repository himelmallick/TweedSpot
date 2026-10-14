# Group-specific (for example cell-type-specific) spatial testing.
#
# With `tweedspot(group = "cell_type")`, each feature is tested for spatial
# variation *within* each group: the locations of one group are fitted with the
# ordinary TweedSpot model, so a feature whose expression differs between
# groups but is flat within each group (a pure composition effect) is not
# called spatially variable. Size factors and the coordinate scaling are
# computed once on all locations, so every group shares one scale.

tweedspot_group_adjust_methods <- function() {
  c("stagewise", "global", "per_group")
}

tweedspot_group_labels <- function(input, group) {
  cd <- SummarizedExperiment::colData(input)
  if (!is.character(group) || length(group) != 1L || !group %in% colnames(cd)) {
    stop("`group` must name a column of `colData(input)`, for example \"cell_type\".")
  }
  labels <- cd[[group]]
  if (!(is.character(labels) || is.factor(labels) || is.logical(labels) ||
        (is.numeric(labels) && all(labels == round(labels), na.rm = TRUE)))) {
    stop("`colData(input)$", group, "` must hold discrete labels (character, factor, ",
         "logical, or integer codes).")
  }
  labels <- droplevels(as.factor(labels))
  if (!nlevels(labels)) {
    stop("`colData(input)$", group, "` has no non-missing labels.")
  }
  labels
}

# Subset the covariate data to the locations of one group. Covariates that are
# constant within the group cannot be estimated there and are dropped, with
# every model term that involves them.
tweedspot_subset_covariates <- function(covariates, idx, group_level, verbose = FALSE) {
  if (is.null(covariates)) {
    return(NULL)
  }
  data <- covariates$data[idx, , drop = FALSE]
  constant <- vapply(data, function(x) length(unique(x)) < 2L, logical(1))
  if (!any(constant)) {
    return(list(data = data, terms = covariates$terms))
  }
  if (verbose) {
    message("Group '", group_level, "': dropping covariate(s) constant within the group: ",
            paste(names(data)[constant], collapse = ", "), ".")
  }
  tt <- stats::terms(stats::as.formula(paste("~", covariates$terms)))
  labels <- attr(tt, "term.labels")
  uses_constant <- vapply(labels, function(lab) {
    any(all.vars(stats::as.formula(paste("~", lab))) %in% names(data)[constant])
  }, logical(1))
  keep_labels <- labels[!uses_constant]
  # offset() terms are not in term.labels; keep them (a constant offset is fine).
  offset_idx <- attr(tt, "offset")
  offsets <- if (is.null(offset_idx)) character(0) else {
    vapply(offset_idx, function(i) {
      paste(deparse(attr(tt, "variables")[[i + 1L]], width.cutoff = 500L), collapse = " ")
    }, character(1))
  }
  keep_labels <- c(keep_labels, offsets)
  if (!length(keep_labels)) {
    return(NULL)
  }
  keep_vars <- unique(unlist(lapply(keep_labels, function(lab) {
    all.vars(stats::as.formula(paste("~", lab)))
  })))
  for (j in seq_len(ncol(data))) {
    if (is.factor(data[[j]])) {
      data[[j]] <- droplevels(data[[j]])
    }
  }
  list(data = data[, names(data) %in% keep_vars, drop = FALSE],
       terms = paste(keep_labels, collapse = " + "))
}

tweedspot_group_fit <- function(fit_input, labels, min_group_size, two_part,
                                family, fit_method, use_bam, bam_discrete,
                                bam_nthreads, smooth_k, BPPARAM, p_criterion,
                                verbose = FALSE) {
  Y <- fit_input$Y
  features <- rownames(Y)
  if (is.null(features)) {
    features <- paste0("feature", seq_len(nrow(Y)))
  }
  per_group <- lapply(levels(labels), function(g) {
    idx <- which(labels == g)
    n_g <- length(idx)
    Yg <- Y[, idx, drop = FALSE]
    n_detected <- rowSums(Yg > 0)
    out <- data.frame(
      feature = features,
      feature_index = seq_len(nrow(Y)),
      group = g,
      n_locations = n_g,
      n_detected = as.integer(n_detected),
      tested = FALSE,
      reason = "",
      stat = NA_real_,
      pval = NA_real_,
      edf = NA_real_,
      dev_expl = NA_real_,
      p = NA_real_,
      stringsAsFactors = FALSE
    )
    if (n_g < min_group_size) {
      out$reason <- "group_too_small"
      return(out)
    }
    coords_g <- fit_input$coords[idx, , drop = FALSE]
    min_support <- tweedspot_basis_k(coords_g, smooth_k)
    supported <- n_detected >= min_support
    out$reason[!supported] <- "low_model_support"
    if (!any(supported)) {
      return(out)
    }
    libsz <- fit_input$libsz
    libsz_g <- if (is.matrix(libsz)) libsz[supported, idx, drop = FALSE] else libsz[idx]
    cov_g <- tweedspot_subset_covariates(fit_input$covariates, idx, g, verbose)
    if (verbose) {
      message("Group '", g, "': testing ", sum(supported), " feature(s) across ",
              n_g, " location(s).")
    }
    res <- tryCatch(
      tweedspot_agnostic(
        Y = Yg[supported, , drop = FALSE],
        coords = coords_g,
        libsz = libsz_g,
        covariates = cov_g,
        two_part = two_part,
        family = family,
        fit_method = fit_method,
        use_bam = use_bam,
        bam_discrete = bam_discrete,
        bam_nthreads = bam_nthreads,
        smooth_k = smooth_k,
        BPPARAM = BPPARAM,
        p_criterion = p_criterion
      ),
      error = function(e) {
        warning("All TweedSpot fits failed in group '", g, "': ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) {
      out$reason[supported] <- "fit_failed"
      return(out)
    }
    out$stat[supported] <- res$stat
    out$pval[supported] <- res$pval
    out$edf[supported] <- res$edf
    out$dev_expl[supported] <- res$dev_expl
    out$p[supported] <- res$p
    fitted_ok <- supported & !is.na(out$pval)
    out$tested[fitted_ok] <- TRUE
    out$reason[supported & !fitted_ok] <- "fit_failed"
    out
  })
  res <- do.call(rbind, per_group)
  res$family <- tweedspot_power_label(res$p)
  rownames(res) <- NULL
  res
}

# Multiple-testing adjustment of feature x group p-values.
#   "stagewise" (Heller et al. 2009; Van den Berge et al. 2017, stageR):
#     1. Screen: combine each feature's group p-values with the Cauchy
#        combination test and adjust across the m screened features with
#        `padj_method`. R features pass at level alpha.
#     2. Confirm: within each feature, Holm-adjust its group p-values at the
#        reduced level alpha * R / m, which accounts for confirming only the
#        features selected in the screen.
#     A pair's reported adjusted p-value is max(screen padj, Holm p * m / R),
#     so a pair is significant at `alpha` exactly when its feature passes the
#     screen and the group passes confirmation. These values are calibrated
#     for the `alpha` used in the run, not for other thresholds.
#   "global": `padj_method` across all tested pairs.
#   "per_group": `padj_method` within each group.
tweedspot_group_adjust <- function(res, method, padj_method, alpha = 0.05) {
  res$padj <- NA_real_
  res$padj_feature <- NA_real_
  ok <- res$tested & !is.na(res$pval)
  if (!any(ok)) {
    return(res)
  }
  if (identical(method, "global")) {
    res$padj[ok] <- stats::p.adjust(res$pval[ok], method = padj_method)
    return(res)
  }
  if (identical(method, "per_group")) {
    for (g in unique(res$group[ok])) {
      sel <- ok & res$group == g
      res$padj[sel] <- stats::p.adjust(res$pval[sel], method = padj_method)
    }
    return(res)
  }
  # Key by row index of the feature, not by name, so duplicated names stay apart.
  key <- res$feature_index[ok]
  rows_by_feature <- split(which(ok), key)
  screen <- vapply(rows_by_feature, function(rows) cct(res$pval[rows]), numeric(1))
  screen_adj <- stats::p.adjust(screen, method = padj_method)
  m <- length(screen_adj)
  n_pass <- sum(screen_adj <= alpha)
  inflate <- m / max(n_pass, 1L)
  for (j in seq_along(rows_by_feature)) {
    rows <- rows_by_feature[[j]]
    holm <- pmin(1, stats::p.adjust(res$pval[rows], method = "holm") * inflate)
    res$padj_feature[rows] <- screen_adj[[j]]
    res$padj[rows] <- pmax(holm, screen_adj[[j]])
  }
  res
}

tweedspot_store_group_results <- function(input, res, group, group_adjust, alpha) {
  name_col <- tweedspot_name_cols(input)
  if (length(name_col)) {
    nm <- as.character(SummarizedExperiment::rowData(input)[[name_col[1L]]])
    res$feature_name <- nm[res$feature_index]
    res <- res[, c("feature", "feature_name",
                   setdiff(names(res), c("feature", "feature_name")))]
  }
  S4Vectors::metadata(input)$tweedspot_group <- list(
    results = res,
    group = group,
    group_adjust = group_adjust,
    alpha = alpha,
    features = rownames(input),
    locations = colnames(input)
  )
  # Per-feature summaries, computed in one pass keyed by row index.
  n <- nrow(input)
  ok <- !is.na(res$padj)
  feature_padj <- rep(NA_real_, n)
  if (any(ok)) {
    summary_p <- if (identical(group_adjust, "stagewise")) res$padj_feature else res$padj
    per_feature <- tapply(summary_p[ok], res$feature_index[ok], min)
    feature_padj[as.integer(names(per_feature))] <- per_feature
  }
  sig <- rep("", n)
  hit <- ok & res$padj <= alpha
  if (any(hit)) {
    groups_hit <- tapply(res$group[hit], res$feature_index[hit], paste, collapse = ";")
    sig[as.integer(names(groups_hit))] <- groups_hit
  }
  rd <- SummarizedExperiment::rowData(input)
  rd$tweedspot_group_padj <- feature_padj
  rd$tweedspot_group_sig <- unname(sig)
  SummarizedExperiment::rowData(input) <- rd
  input
}
