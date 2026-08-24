#' Detect spatially variable genes with self-adaptive spatial Tweedie models
#'
#' `tweedspot()` is the single entry point for TweedSpot. It dispatches on its
#' arguments to one of three modes:
#' \itemize{
#'   \item **agnostic** (neither `W` nor `celltype` supplied): fits a Tweedie GAM
#'     `s(x, y)` per gene and tests overall spatial variation.
#'   \item **proportion-based, cell-type-specific** (`W` supplied): for spot-level
#'     data, fits a separate spatial smooth per cell type weighted by the
#'     cell-type proportion matrix `W` (in the spirit of CTSV / C-SIDE), and
#'     returns one spatial p-value per cell type.
#'   \item **label-based, cell-type-specific** (`celltype` supplied): for
#'     single-cell-resolution data, fits `s(x, y, by = celltype)` and returns
#'     one spatial p-value per cell type.
#' }
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a count assay and
#'   spatial coordinates.
#' @param assay_name Name of the count assay to model. Default `"counts"`.
#' @param W Optional spots-by-K matrix of cell-type proportions (rows sum to 1).
#'   Triggers proportion-based cell-type-specific detection. Pass true
#'   proportions in simulation; pass deconvolution estimates (e.g. RCTD) on real
#'   data.
#' @param celltype Optional factor of length `ncol(input)` giving a cell-type
#'   label per location. Triggers label-based cell-type-specific detection.
#'   Ignored if `W` is supplied.
#' @param two_part Logical; if `TRUE`, fit the two-part (hurdle) model (Tweedie
#'   on counts plus a binomial model on presence/absence) and combine the two
#'   p-values. Agnostic mode only. Default `FALSE`.
#' @param combine P-value combination rule for the two-part model: `"CCT"`
#'   (Cauchy), `"stouffer"`, or `"min"`. Default `"CCT"`.
#' @param family `mgcv` family for the count model. Default `"tw"` (Tweedie).
#' @param fit_method Smoothing parameter estimation method passed to `mgcv`,
#'   e.g. `"REML"` (default) or `"GCV.Cp"`.
#' @param use_bam Logical; use [mgcv::bam()] instead of [mgcv::gam()] for
#'   scalability on large datasets. Default `FALSE`.
#' @param bam_discrete Logical; when `use_bam = TRUE`, use discretized fitting
#'   via `mgcv::bam(discrete = TRUE)` for speed. Default `TRUE`.
#' @param bam_nthreads Integer; number of threads for `mgcv::bam()`. Default 1.
#' @param smooth_k Optional basis dimension for the spatial smooth. Smaller
#'   values can speed fitting at the cost of flexibility. Default `NULL`, which
#'   uses `mgcv`'s default.
#' @param padj_method Multiple-testing correction passed to [stats::p.adjust()].
#'   Default `"BY"`.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallelization across genes. Default [BiocParallel::bpparam()], which
#'   uses the currently registered BiocParallel backend.
#'
#' @return The input `SpatialExperiment` with results added. In agnostic mode the
#'   per-gene p-value, adjusted p-value, and statistic are written to
#'   [SummarizedExperiment::rowData()]. In cell-type-specific modes a
#'   genes-by-celltypes p-value matrix and its adjusted counterpart are stored in
#'   `rowData()` as nested matrices, and a tidy long `data.frame` is attached to
#'   `metadata(input)$tweedspot`.
#'
#' @examples
#' \dontrun{
#' library(STexampleData)
#' spe <- ST_mouseOB()
#' # agnostic
#' spe <- tweedspot(spe)
#' head(SummarizedExperiment::rowData(spe))
#' # cell-type-specific with proportions
#' spe <- tweedspot(spe, W = my_proportions)
#' }
#' @export
tweedspot <- function(input,
                      assay_name = "counts",
                      W = NULL,
                      celltype = NULL,
                      two_part = FALSE,
                      combine = c("CCT", "stouffer", "min"),
                      family = "tw",
                      fit_method = "REML",
                      use_bam = FALSE,
                      bam_discrete = TRUE,
                      bam_nthreads = 1L,
                      smooth_k = NULL,
                      padj_method = "BY",
                      BPPARAM = BiocParallel::bpparam()) {

  # Restrict the p-value combination rule to the supported options.
  combine <- match.arg(combine)

  # The package is written around SpatialExperiment inputs.
  stopifnot(methods::is(input, "SpatialExperiment"))

  # Pull out the count matrix, standardized spatial coordinates, gene names,
  # and per-location size factors used as a log-offset in the GAMs.
  Y      <- as.matrix(SummarizedExperiment::assay(input, assay_name))  # genes x spots
  coords <- scale(SpatialExperiment::spatialCoords(input))             # mean 0, sd 1
  genes  <- rownames(input)
  libsz  <- .tweedspot_libsize(input, Y)

  # The spatial smooths require at least two coordinate columns.
  if (ncol(coords) < 2) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  if (!is.null(W)) {
    # Coerce to matrix so the downstream algebra and dimension checks behave
    # predictably even if the user supplies a data.frame.
    W <- as.matrix(W)
    if (!is.numeric(W)) {
      stop("`W` must be a numeric spots-by-celltypes matrix.")
    }
    # The proportion matrix is spot-by-celltype, so its row count must match
    # the number of spatial locations in the input object.
    if (nrow(W) != ncol(input)) {
      stop("`W` must have one row per spatial location in `input`.")
    }
    if (ncol(W) < 1L) {
      stop("`W` must contain at least one cell type column.")
    }
    if (anyNA(W) || any(!is.finite(W))) {
      stop("`W` must not contain missing or non-finite values.")
    }
    if (any(W < 0)) {
      stop("`W` must contain non-negative cell-type proportions.")
    }
    if (any(rowSums(W) <= 0)) {
      stop("Each row of `W` must have a positive sum.")
    }
  }
  if (!is.null(celltype)) {
    # Label-based mode expects one cell-type label per spatial location.
    if (length(celltype) != ncol(input)) {
      stop("`celltype` must have length equal to the number of spatial locations in `input`.")
    }
    if (anyNA(celltype)) {
      stop("`celltype` must not contain missing values.")
    }
  }

  # Select the analysis mode from the supplied cell-type information.
  mode <- if (!is.null(W)) "proportion" else if (!is.null(celltype)) "label" else "agnostic"

  # The two-part model is only implemented for the agnostic workflow.
  if (two_part && mode != "agnostic") {
    stop("`two_part = TRUE` is only supported in agnostic mode.")
  }

  if (mode == "agnostic") {
    # Fit one spatial model per gene and return one statistic/p-value per gene.
    res <- .tweedspot_agnostic(Y, coords, libsz, two_part, combine,
                               family, fit_method, use_bam, bam_discrete,
                               bam_nthreads, smooth_k, BPPARAM)

    # Write agnostic-mode results directly into rowData for easy ranking.
    SummarizedExperiment::rowData(input)$tweedspot_stat <- res$stat
    SummarizedExperiment::rowData(input)$tweedspot_pval <- res$pval
    SummarizedExperiment::rowData(input)$tweedspot_padj <-
      stats::p.adjust(res$pval, method = padj_method)
    return(input)
  }

  # Cell-type-specific modes return a genes-by-celltype p-value matrix.
  pmat <- if (mode == "proportion") {
    .tweedspot_proportion(Y, coords, libsz, W, family, fit_method, use_bam,
                          bam_discrete, bam_nthreads, smooth_k, BPPARAM)
  } else {
    .tweedspot_label(Y, coords, libsz, celltype, family, fit_method, use_bam,
                     bam_discrete, bam_nthreads, smooth_k, BPPARAM)
  }

  # Preserve gene IDs on the p-value matrix rows.
  rownames(pmat) <- genes

  # Adjust p-values separately within each cell type.
  padj <- apply(pmat, 2, stats::p.adjust, method = padj_method)

  # Store the wide matrix outputs in rowData.
  rd <- SummarizedExperiment::rowData(input)
  rd$tweedspot_pval <- pmat
  rd$tweedspot_padj <- padj
  SummarizedExperiment::rowData(input) <- rd

  # Also attach a long-format table for downstream plotting/reporting.
  long <- data.frame(
    gene     = rep(genes, times = ncol(pmat)),
    celltype = rep(colnames(pmat), each = nrow(pmat)),
    pval     = as.vector(pmat),
    padj     = as.vector(padj),
    stringsAsFactors = FALSE)
  S4Vectors::metadata(input)$tweedspot <- long
  input
}

# ---- internal: library size (scran if available, else total-count) ----------
.tweedspot_libsize <- function(input, Y) {
  # Reuse existing size factors if the input object already has them, since that
  # avoids recomputation and lets upstream preprocessing control normalization.
  if (methods::is(input, "SingleCellExperiment")) {
    sf <- SingleCellExperiment::sizeFactors(input)
    if (!is.null(sf) && length(sf) == ncol(input) && all(is.finite(sf) | is.na(sf))) {
      sf[is.na(sf) | sf <= 0] <- 1
      return(sf)
    }
  }
  if (requireNamespace("scran", quietly = TRUE) &&
      requireNamespace("scuttle", quietly = TRUE)) {
    # Prefer deconvolution-based size factors when the Bioconductor helpers
    # are installed; otherwise fall back to simple total-count normalization.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Y))
    sce <- scran::computeSumFactors(sce, BPPARAM = BiocParallel::SerialParam())
    sf  <- SingleCellExperiment::sizeFactors(sce); sf[is.na(sf) | sf <= 0] <- 1
    return(sf)
  }

  # Last-resort fallback: normalize by total counts relative to the median spot.
  ls <- colSums(Y); ls / stats::median(ls)
}

.tweedspot_fit <- function(formula, family, data, fit_method, use_bam,
                           bam_discrete, bam_nthreads) {
  # mgcv's Tweedie family behaves more reliably when the namespace is attached,
  # so do that lazily only for the Tweedie path.
  if (identical(family, "tw") && !"package:mgcv" %in% search()) {
    base::attachNamespace(asNamespace("mgcv"))
  }

  # Use bam() when requested for larger datasets; otherwise use gam().
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

.tweedspot_basis_k <- function(coords, smooth_k = NULL, default_k = 30L) {
  # mgcv smooths cannot use more basis functions than the spatial resolution can
  # support, so cap k by the number of unique coordinate pairs.
  n_unique <- nrow(unique(as.data.frame(coords[, seq_len(min(2, ncol(coords))), drop = FALSE])))
  if (n_unique <= 3L) {
    return(3L)
  }
  if (is.null(smooth_k)) {
    return(as.integer(min(default_k, n_unique - 1L)))
  }
  as.integer(max(3L, min(as.integer(smooth_k), n_unique - 1L)))
}

# ---- internal: agnostic per-gene fit ----------------------------------------
.tweedspot_agnostic <- function(Y, coords, libsz, two_part, combine,
                                family, fit_method, use_bam, bam_discrete,
                                bam_nthreads, smooth_k, BPPARAM) {
  # Base design data shared by every gene-specific fit.
  d_base <- data.frame(x = coords[, 1], y = coords[, 2], libsz = libsz)

  # Choose a stable basis dimension for the spatial smooth.
  smooth_k <- .tweedspot_basis_k(coords, smooth_k)

  # Build the main count-model formula and the optional presence/absence formula.
  smooth_term <- sprintf("s(x, y, k = %d)", smooth_k)
  form <- stats::as.formula(sprintf("expr ~ %s + offset(log(libsz))", smooth_term))
  form_bin <- stats::as.formula(sprintf("I(expr > 0) ~ %s + offset(log(libsz))", smooth_term))

  one <- function(i) {
    # Pull one gene across all spatial locations into the modeling frame.
    d <- cbind(expr = as.numeric(Y[i, ]), d_base)

    # Fit one spatial smooth per gene and test whether the smooth term is
    # needed after accounting for library size.
    f <- tryCatch(suppressWarnings(
      .tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete,
                     bam_nthreads)), error = function(e) NULL)

    # Failed fits are returned as missing results so the rest of the run can continue.
    if (is.null(f)) return(c(stat = NA_real_, pval = NA_real_, plogit = NA_real_))

    # mgcv stores smooth-term inference in s.table; use the first smooth row
    # because there is only one spatial smooth in agnostic mode.
    s <- summary(f)$s.table
    stat <- s[1, 3]; pv <- s[1, "p-value"]; plog <- NA_real_
    if (two_part) {
      # Optional hurdle-style companion model on detection/non-detection.
      f1 <- tryCatch(suppressWarnings(
        .tweedspot_fit(form_bin, "binomial", d, fit_method, use_bam,
                       bam_discrete, bam_nthreads)), error = function(e) NULL)

      # If the binomial fit succeeds, pull its smooth-term p-value as well.
      if (!is.null(f1)) plog <- summary(f1)$s.table[1, "p-value"]
    }
    c(stat = stat, pval = pv, plogit = plog)
  }

  # Parallelize across genes using the user-supplied BiocParallel backend.
  M <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(Y)), one, BPPARAM = BPPARAM))

  # In the one-part model the count-model p-value is the final p-value.
  if (!two_part) return(list(stat = M[, "stat"], pval = M[, "pval"]))

  # In the two-part model combine the count and detection p-values per gene.
  combined <- switch(combine,
    CCT      = vapply(seq_len(nrow(M)), function(i)
                 .CCT(c(M[i, "pval"], M[i, "plogit"])), numeric(1)),
    stouffer = vapply(seq_len(nrow(M)), function(i) {
                 p <- c(M[i, "pval"], M[i, "plogit"]); p <- p[!is.na(p)]
                 if (!length(p)) return(NA_real_)
                 1 - stats::pnorm(sum(stats::qnorm(1 - p)) / sqrt(length(p))) }, numeric(1)),
    min      = pmin(M[, "pval"], M[, "plogit"], na.rm = TRUE))
  list(stat = M[, "stat"], pval = combined)
}

# ---- internal: proportion-based (CTSV-style) --------------------------------
.tweedspot_proportion <- function(Y, coords, libsz, W, family, fit_method, use_bam,
                                  bam_discrete, bam_nthreads, smooth_k, BPPARAM) {
  # Normalize rows so each spot's proportions sum to one.
  K  <- ncol(W); W <- W / rowSums(W)

  # Pull x/y coordinates once and reuse them across genes.
  h1 <- coords[, 1]; h2 <- coords[, 2]
  smooth_k <- .tweedspot_basis_k(coords, smooth_k)
  # Build cell-type-weighted coordinates so each smooth represents the spatial
  # pattern attributable to one cell type in mixed spots.
  Xx <- W * h1; Xy <- W * h2
  colnames(Xx) <- paste0("ct", seq_len(K), "_x")
  colnames(Xy) <- paste0("ct", seq_len(K), "_y")
  colnames(W)  <- if (is.null(colnames(W))) paste0("ct", seq_len(K)) else colnames(W)
  smooth_spec <- function(xcol, ycol) {
    # Each cell type gets its own spatial smooth over weighted coordinates.
    sprintf("s(%s, %s, k = %d)", xcol, ycol, smooth_k)
  }

  # The model includes one smooth per cell type plus linear terms for the
  # proportions themselves, with no intercept.
  form <- stats::as.formula(sprintf(
    "z ~ -1 + offset(log(libsz)) + %s + %s",
    paste(mapply(smooth_spec, colnames(Xx), colnames(Xy), USE.NAMES = FALSE), collapse = " + "),
    paste(colnames(W), collapse = " + ")))

  one <- function(i) {
    # Bind the gene's counts with the precomputed design matrices.
    dat <- cbind.data.frame(z = as.numeric(Y[i, ]), Xx, Xy, W, libsz = libsz)
    # Each smooth term corresponds to one cell type and yields one p-value.
    f <- tryCatch(suppressWarnings(
      .tweedspot_fit(form, family, dat, fit_method, use_bam, bam_discrete,
                     bam_nthreads)), error = function(e) NULL)
    if (is.null(f)) return(rep(NA_real_, K))

    # summary(f)$s.table has one row per smooth; align those rows to the
    # cell-type columns and leave anything unmatched as NA.
    st <- summary(f)$s.table; out <- rep(NA_real_, K)
    out[seq_len(min(K, nrow(st)))] <- st[seq_len(min(K, nrow(st))), "p-value"]
    out
  }

  # Parallelize the cell-type-specific fits across genes.
  M <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(Y)), one, BPPARAM = BPPARAM))
  colnames(M) <- colnames(W); M
}

# ---- internal: label-based (single-cell resolution) -------------------------
.tweedspot_label <- function(Y, coords, libsz, celltype, family, fit_method, use_bam,
                             bam_discrete, bam_nthreads, smooth_k, BPPARAM) {
  # Factor encoding defines the set of cell-type-specific smooths to fit.
  celltype <- as.factor(celltype)
  lvls     <- levels(celltype)

  # Base design shared across all per-gene fits in label-based mode.
  d_base   <- data.frame(x = coords[, 1], y = coords[, 2], celltype = celltype, libsz = libsz)
  smooth_k <- .tweedspot_basis_k(coords, smooth_k)

  # by = celltype asks mgcv for a separate spatial smooth for each label.
  smooth_term <- sprintf("s(x, y, by = celltype, k = %d)", smooth_k)
  form <- stats::as.formula(sprintf("expr ~ %s + offset(log(libsz))", smooth_term))

  one <- function(i) {
    # Pull one gene's counts into the shared design frame.
    d <- cbind(expr = as.numeric(Y[i, ]), d_base)
    # Fit separate spatial smooths for each cell-type label in a
    # single-cell-resolution setting.
    f <- tryCatch(suppressWarnings(
      .tweedspot_fit(form, family, d, fit_method, use_bam, bam_discrete,
                     bam_nthreads)), error = function(e) NULL)
    if (is.null(f)) return(stats::setNames(rep(NA_real_, length(lvls)), lvls))

    # Map the returned smooth-term p-values back to the factor levels.
    st <- summary(f)$s.table
    p  <- stats::setNames(rep(NA_real_, length(lvls)), lvls)
    p[seq_len(min(length(lvls), nrow(st)))] <- st[seq_len(min(length(lvls), nrow(st))), "p-value"]
    p
  }

  # Parallelize across genes and return a genes-by-celltype p-value matrix.
  M <- do.call(rbind, BiocParallel::bplapply(seq_len(nrow(Y)), one, BPPARAM = BPPARAM))
  colnames(M) <- lvls; M
}

# ---- internal: Cauchy combination (carried over from recovered code) --------
.CCT <- function(pvals, weights = NULL) {
  # Replace NA with 1 so missing components contribute no signal.
  pvals <- ifelse(is.na(pvals), 1, pvals)

  # Clamp values away from exactly 0 or 1 for numerical stability.
  pvals <- pmin(pmax(pvals, .Machine$double.xmin), 1 - .Machine$double.eps)

  # Default to equal weights, or renormalize any user-supplied weights.
  if (is.null(weights)) weights <- rep(1 / length(pvals), length(pvals))
  else weights <- weights / sum(weights)

  # Extremely tiny p-values are handled with the asymptotic form to avoid tan()
  # overflow; otherwise use the standard Cauchy transform.
  is.small <- pvals < 1e-16
  cct <- if (!any(is.small)) sum(weights * tan((0.5 - pvals) * pi))
         else sum((weights[is.small] / pvals[is.small]) / pi) +
              sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))

  # Convert the Cauchy statistic back to a p-value.
  if (cct > 1e15) (1 / cct) / pi else 1 - stats::pcauchy(cct)
}
