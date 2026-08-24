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

tweedspot_covariates <- function(input, covariates) {
  if (is.null(covariates)) {
    return(NULL)
  }

  coldata <- as.data.frame(SummarizedExperiment::colData(input),
                           stringsAsFactors = FALSE)

  if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L) {
      stop("`covariates` must be a one-sided formula like `~ batch + age`.")
    }
    covariate_terms <- paste(deparse(covariates[[2L]]), collapse = "")
    covariate_vars <- all.vars(covariates)
  } else if (is.character(covariates)) {
    covariate_vars <- covariates
    covariate_terms <- paste(covariate_vars, collapse = " + ")
  } else {
    stop("`covariates` must be NULL, a one-sided formula, or a character vector of `colData(input)` column names.")
  }

  if (!length(covariate_vars)) {
    stop("`covariates` must reference at least one column in `colData(input)`.")
  }

  missing_vars <- setdiff(covariate_vars, colnames(coldata))
  if (length(missing_vars)) {
    stop(sprintf("`covariates` references columns not found in `colData(input)`: %s.",
                 paste(missing_vars, collapse = ", ")))
  }

  covariate_data <- coldata[, covariate_vars, drop = FALSE]
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
  if (any(!is.finite(libsz) | libsz <= 0)) {
    stop(
      "Library sizes must be positive and finite for every spatial location. ",
      "This often happens when filtered spots have zero total counts; remove those spots ",
      "or provide a positive `colData(input)$sizeFactor`."
    )
  }
  d_base <- data.frame(x = coords[, 1], y = coords[, 2], libsz = libsz)
  if (!is.null(covariates)) {
    d_base <- cbind(d_base, covariates$data)
  }
  smooth_k <- tweedspot_basis_k(coords, smooth_k)
  smooth_term <- sprintf("s(x, y, k = %d)", smooth_k)
  cov_terms <- if (is.null(covariates)) character(0) else covariates$terms
  rhs <- paste(c(smooth_term, cov_terms, "offset(log(libsz))"), collapse = " + ")
  form <- stats::as.formula(sprintf("expr ~ %s", rhs))
  d <- cbind(expr = as.numeric(expr), d_base)
  tryCatch(
    suppressWarnings(
      tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete, bam_nthreads)
    ),
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

tweedspot_agnostic <- function(Y, coords, libsz, covariates, two_part, combine,
                               family, fit_method, use_bam, bam_discrete,
                               bam_nthreads, smooth_k, BPPARAM) {
  if (any(!is.finite(libsz) | libsz <= 0)) {
    stop(
      "Library sizes must be positive and finite for every spatial location. ",
      "This often happens when filtered spots have zero total counts; remove those spots ",
      "or provide a positive `colData(input)$sizeFactor`."
    )
  }
  d_base <- data.frame(x = coords[, 1], y = coords[, 2], libsz = libsz)
  if (!is.null(covariates)) {
    d_base <- cbind(d_base, covariates$data)
  }
  smooth_k <- tweedspot_basis_k(coords, smooth_k)

  smooth_term <- sprintf("s(x, y, k = %d)", smooth_k)
  cov_terms <- if (is.null(covariates)) character(0) else covariates$terms
  rhs <- paste(c(smooth_term, cov_terms, "offset(log(libsz))"), collapse = " + ")
  form <- stats::as.formula(sprintf("expr ~ %s", rhs))
  form_bin <- stats::as.formula(sprintf("I(expr > 0) ~ %s", rhs))

  one <- function(i) {
    d <- cbind(expr = as.numeric(Y[i, ]), d_base)
    f <- tryCatch(suppressWarnings(
      tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete,
                    bam_nthreads)), error = function(e) NULL)

    if (is.null(f)) {
      return(c(
        stat = NA_real_,
        pval = NA_real_,
        plogit = NA_real_,
        edf = NA_real_,
        dev_expl = NA_real_
      ))
    }

    s <- summary(f)$s.table
    model_summary <- summary(f)
    stat <- s[1, 3]
    pv <- s[1, "p-value"]
    plog <- NA_real_
    edf <- s[1, 1]
    dev_expl <- model_summary$dev.expl
    if (two_part) {
      f1 <- tryCatch(suppressWarnings(
        tweedspot_fit(form_bin, "binomial", d, fit_method, use_bam,
                      bam_discrete, bam_nthreads)), error = function(e) NULL)

      if (!is.null(f1)) plog <- summary(f1)$s.table[1, "p-value"]
    }
    c(stat = stat, pval = pv, plogit = plog, edf = edf, dev_expl = dev_expl)
  }

  M <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(Y)), one, BPPARAM = BPPARAM))

  if (!two_part) {
    if (all(is.na(M[, "pval"]))) {
      stop(
        "All gene-wise TweedSpot fits failed. Check covariates, library sizes, and filtering thresholds."
      )
    }
    return(list(
      stat = M[, "stat"],
      pval = M[, "pval"],
      edf = M[, "edf"],
      dev_expl = M[, "dev_expl"]
    ))
  }

  combined <- switch(combine,
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
