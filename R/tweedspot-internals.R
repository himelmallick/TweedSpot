tweedspot_libsize <- function(input, Y) {
  coldata <- SummarizedExperiment::colData(input)
  if ("sizeFactor" %in% names(coldata)) {
    sf <- as.numeric(coldata[["sizeFactor"]])
    if (length(sf) == ncol(input) && all(is.finite(sf) | is.na(sf))) {
      sf[is.na(sf) | sf <= 0] <- 1
      return(sf)
    }
  }

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

tweedspot_resolve_gene <- function(input, gene) {
  if (is.numeric(gene) && length(gene) == 1L && !is.na(gene)) {
    idx <- as.integer(gene)
    if (idx < 1L || idx > nrow(input)) {
      stop("Numeric `gene` index is out of bounds.")
    }
  } else if (is.character(gene) && length(gene) == 1L && !is.na(gene)) {
    idx <- match(gene, rownames(input))
    if (is.na(idx) && "gene_name" %in% colnames(SummarizedExperiment::rowData(input))) {
      idx <- match(gene, as.character(SummarizedExperiment::rowData(input)$gene_name))
    }
    if (is.na(idx)) {
      stop("`gene` was not found among row names or `rowData(input)$gene_name`.")
    }
  } else {
    stop("`gene` must be a single row index, row name, or gene name.")
  }

  label <- rownames(input)[idx]
  if ("gene_name" %in% colnames(SummarizedExperiment::rowData(input))) {
    gene_name <- as.character(SummarizedExperiment::rowData(input)$gene_name[idx])
    if (!is.na(gene_name) && nzchar(gene_name)) {
      label <- paste0(gene_name, " (", rownames(input)[idx], ")")
    }
  }
  list(index = idx, label = label)
}

tweedspot_plot_coords <- function(input) {
  coords <- SpatialExperiment::spatialCoords(input)
  if (ncol(coords) < 2L) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  coords[, seq_len(2), drop = FALSE]
}

tweedspot_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#D9D4CC", fill = NA, linewidth = 0.6),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

tweedspot_annotation_exclusions <- function(input, exclude_mito, exclude_ribo,
                                            gene_name_col = NULL,
                                            gene_biotype_col = NULL) {
  rd <- as.data.frame(SummarizedExperiment::rowData(input), stringsAsFactors = FALSE)
  gene_names <- rownames(input)
  biotypes <- NULL

  if (!is.null(gene_name_col)) {
    if (!gene_name_col %in% colnames(rd)) {
      stop(sprintf("`gene_name_col` not found in `rowData(input)`: %s.", gene_name_col))
    }
    gene_names <- as.character(rd[[gene_name_col]])
  }
  if (!is.null(gene_biotype_col)) {
    if (!gene_biotype_col %in% colnames(rd)) {
      stop(sprintf("`gene_biotype_col` not found in `rowData(input)`: %s.", gene_biotype_col))
    }
    biotypes <- as.character(rd[[gene_biotype_col]])
  }

  gene_names[is.na(gene_names)] <- ""
  mito <- rep(FALSE, nrow(input))
  ribo <- rep(FALSE, nrow(input))

  if (exclude_mito) {
    mito_name <- grepl("^mt[-.]", gene_names, ignore.case = TRUE) |
      grepl("^mt[a-z0-9]", gene_names, ignore.case = TRUE)
    mito_biotype <- if (is.null(biotypes)) rep(FALSE, nrow(input)) else {
      grepl("mitochond", biotypes, ignore.case = TRUE)
    }
    mito <- mito_name | mito_biotype
  }

  if (exclude_ribo) {
    ribo_name <- grepl("^(rpl|rps|mrpl|mrps)", gene_names, ignore.case = TRUE)
    ribo_biotype <- if (is.null(biotypes)) rep(FALSE, nrow(input)) else {
      grepl("ribosom", biotypes, ignore.case = TRUE)
    }
    ribo <- ribo_name | ribo_biotype
  }

  mito | ribo
}

tweedspot_gene_filters <- function(Y, filter_genes_ncounts, filter_genes_pcspots,
                                   filter_genes_nspots, filter_genes_mean,
                                   filter_genes_var) {
  keep_ncounts <- rowSums(Y) >= filter_genes_ncounts
  keep_pcspots <- rowMeans(Y > 0) * 100 >= filter_genes_pcspots
  keep_nspots <- rowSums(Y > 0) >= filter_genes_nspots
  keep_mean <- rowMeans(Y) >= filter_genes_mean
  keep_var <- apply(Y, 1, stats::var) >= filter_genes_var
  list(
    keep = keep_ncounts & keep_pcspots & keep_nspots & keep_mean & keep_var,
    keep_ncounts = keep_ncounts,
    keep_pcspots = keep_pcspots,
    keep_nspots = keep_nspots,
    keep_mean = keep_mean,
    keep_var = keep_var
  )
}

tweedspot_report_filtering <- function(output, n_genes_before, n_spots_before,
                                       criteria_removed, dropped_zero_spots) {
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

tweedspot_parse_covariates <- function(covariates) {
  if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L) {
      stop("`covariates` must be a one-sided formula like `~ batch + age`.")
    }
    return(list(
      vars = all.vars(covariates),
      terms = paste(deparse(covariates[[2L]]), collapse = "")
    ))
  }
  if (is.character(covariates)) {
    return(list(
      vars = covariates,
      terms = paste(covariates, collapse = " + ")
    ))
  }
  stop("`covariates` must be NULL, a one-sided formula, or a character vector of `colData(input)` column names.")
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

tweedspot_sanitize_covariates <- function(covariate_data, covariate_vars, covariate_terms) {
  reserved <- c("expr", "x", "y", "libsz")
  sanitized_names <- make.names(names(covariate_data), unique = TRUE)
  if (any(sanitized_names %in% reserved)) {
    stop("Referenced `covariates` columns must not use reserved names: expr, x, y, libsz.")
  }
  names(covariate_data) <- sanitized_names
  for (j in seq_along(covariate_vars)) {
    covariate_terms <- gsub(
      paste0("\\b", covariate_vars[j], "\\b"),
      sanitized_names[j],
      covariate_terms
    )
  }
  list(data = covariate_data, terms = covariate_terms)
}

tweedspot_covariates <- function(input, covariates) {
  if (is.null(covariates)) {
    return(NULL)
  }

  coldata <- as.data.frame(SummarizedExperiment::colData(input),
                           stringsAsFactors = FALSE)
  parsed <- tweedspot_parse_covariates(covariates)
  covariate_vars <- parsed$vars
  covariate_terms <- parsed$terms
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
  tweedspot_sanitize_covariates(covariate_data, covariate_vars, covariate_terms)
}

tweedspot_fit <- function(formula, family, data, fit_method, use_bam,
                          bam_discrete, bam_nthreads) {
  if (identical(family, "tw") && !"package:mgcv" %in% search()) {
    base::attachNamespace(asNamespace("mgcv"))
  }

  if (use_bam) {
    return(mgcv::bam(
      formula = formula,
      family = family,
      data = data,
      method = fit_method,
      discrete = bam_discrete,
      nthreads = bam_nthreads
    ))
  }
  mgcv::gam(
    formula = formula,
    family = family,
    data = data,
    method = fit_method
  )
}

tweedspot_fit_single_gene <- function(expr, coords, libsz, covariates, family,
                                      fit_method, use_bam, bam_discrete,
                                      bam_nthreads, smooth_k) {
  tweedspot_validate_libsize(libsz)
  formula_parts <- tweedspot_model_formulas(coords, libsz, covariates, smooth_k)
  form <- formula_parts$form
  d_base <- formula_parts$data
  d <- cbind(expr = as.numeric(expr), d_base)
  tryCatch(
    tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete, bam_nthreads),
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
  d_base <- data.frame(x = coords[, 1], y = coords[, 2], libsz = libsz)
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

tweedspot_smooth_summary <- function(fit, two_part, form_bin, data, fit_method,
                                     use_bam, bam_discrete, bam_nthreads) {
  s <- summary(fit)$s.table
  model_summary <- summary(fit)
  plog <- NA_real_
  if (two_part) {
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
    dev_expl = model_summary$dev.expl
  )
}

tweedspot_fit_row <- function(i, Y, d_base, form, form_bin, two_part, family,
                              fit_method, use_bam, bam_discrete, bam_nthreads) {
  d <- cbind(expr = as.numeric(Y[i, ]), d_base)
  fit <- tryCatch(
    tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete, bam_nthreads),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(c(
      stat = NA_real_,
      pval = NA_real_,
      plogit = NA_real_,
      edf = NA_real_,
      dev_expl = NA_real_
    ))
  }
  tweedspot_smooth_summary(
    fit = fit,
    two_part = two_part,
    form_bin = form_bin,
    data = d,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    bam_nthreads = bam_nthreads
  )
}

tweedspot_collect_fits <- function(Y, d_base, form, form_bin, two_part, family,
                                   fit_method, use_bam, bam_discrete,
                                   bam_nthreads, BPPARAM) {
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
      bam_nthreads = bam_nthreads
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
    dev_expl = M[, "dev_expl"]
  )
}

tweedspot_combine_pvalues <- function(M, combine) {
  switch(combine,
    CCT = vapply(seq_len(nrow(M)), function(i) {
      cct(c(M[i, "pval"], M[i, "plogit"]))
    }, numeric(1)),
    stouffer = vapply(seq_len(nrow(M)), function(i) {
      p <- c(M[i, "pval"], M[i, "plogit"])
      p <- p[!is.na(p)]
      if (!length(p)) return(NA_real_)
      1 - stats::pnorm(sum(stats::qnorm(1 - p)) / sqrt(length(p)))
    }, numeric(1)),
    min = pmin(M[, "pval"], M[, "plogit"], na.rm = TRUE)
  )
}

tweedspot_extract_two_part <- function(M, combine) {
  combined <- tweedspot_combine_pvalues(M, combine)
  if (all(is.na(combined))) {
    stop(
      "All gene-wise TweedSpot fits failed. Check covariates, library sizes, and filtering thresholds."
    )
  }
  list(
    stat = M[, "stat"],
    pval = combined,
    edf = M[, "edf"],
    dev_expl = M[, "dev_expl"]
  )
}

tweedspot_agnostic <- function(Y, coords, libsz, covariates, two_part, combine,
                               family, fit_method, use_bam, bam_discrete,
                               bam_nthreads, smooth_k, BPPARAM) {
  tweedspot_validate_libsize(libsz)
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
    BPPARAM = BPPARAM
  )
  if (!two_part) {
    return(tweedspot_extract_one_part(M))
  }
  tweedspot_extract_two_part(M, combine)
}

tweedspot_prepare_input <- function(input, assay_name, covariates) {
  Y <- as.matrix(SummarizedExperiment::assay(input, assay_name))
  coords <- scale(SpatialExperiment::spatialCoords(input))
  if (ncol(coords) < 2L) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  list(
    Y = Y,
    coords = coords,
    libsz = tweedspot_libsize(input, Y),
    covariates = tweedspot_covariates(input, covariates)
  )
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
  input
}

tweedspot_message_end <- function(verbose, pval) {
  if (verbose) {
    message("Completed TweedSpot fits for ", sum(!is.na(pval)), " gene(s).")
  }
}

tweedspot_prepare_plot_fit_data <- function(input, gene, assay_name, covariates,
                                            family, fit_method, use_bam,
                                            bam_discrete, bam_nthreads,
                                            smooth_k) {
  gene_info <- tweedspot_resolve_gene(input, gene)
  Y <- as.matrix(SummarizedExperiment::assay(input, assay_name))
  fit <- tweedspot_fit_single_gene(
    expr = as.numeric(Y[gene_info$index, ]),
    coords = scale(SpatialExperiment::spatialCoords(input)),
    libsz = tweedspot_libsize(input, Y),
    covariates = tweedspot_covariates(input, covariates),
    family = family,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    bam_nthreads = bam_nthreads,
    smooth_k = smooth_k
  )
  if (is.null(fit)) {
    stop("The single-gene spatial model failed to fit for the requested gene.")
  }
  raw_coords <- tweedspot_plot_coords(input)
  list(
    label = gene_info$label,
    df = data.frame(
      x = raw_coords[, 1],
      y = raw_coords[, 2],
      value = as.numeric(stats::predict(fit, type = "terms")[, 1])
    )
  )
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
