#' Detect spatially variable genes with self-adaptive spatial Tweedie models
#'
#' `tweedspot()` fits a Tweedie GAM with a spatial smooth `s(x, y)` per gene
#' and tests each gene for overall spatial variation.
#'
#' @param input A [SpatialExperiment::SpatialExperiment] with a count assay and
#'   spatial coordinates.
#' @param assay_name Name of the count assay to model. Default `"counts"`.
#' @param covariates Optional one-sided formula specifying columns from
#'   `colData(input)` to adjust for, for example `~ batch + age + sex`.
#'   Character vectors of `colData(input)` column names are also supported.
#' @param two_part Logical; if `TRUE`, fit the two-part (hurdle) model (Tweedie
#'   on counts plus a binomial model on presence/absence) and combine the two
#'   p-values with the Cauchy combination test (CCT), which tests whether
#'   either part varies in space and stays valid when the two p-values are
#'   dependent. When the presence/absence part is skipped, the count part's
#'   p-value is used. Default `FALSE`.
#' @param family `mgcv` family for the count model. Default `"tw"` (Tweedie).
#'   Ignored when `tweedie_p` is supplied.
#' @param tweedie_p Tweedie variance power \eqn{p}, where
#'   \eqn{Var(Y) = \phi \mu^p}. If `NULL` (default), TweedSpot estimates
#'   \eqn{p \in (1, 2)} per gene with [mgcv::tw()]. Supported fixed values are
#'   `0` for Gaussian, `1` for Poisson, values in `(1, 2]` for fixed Tweedie
#'   (`2` is Gamma and requires strictly positive responses), and `3` for
#'   inverse Gaussian with strictly positive responses. All families use a log
#'   link so the size-factor offset acts multiplicatively. Use `"auto"` to
#'   choose the variance power per gene by an information criterion (see
#'   `p_criterion` and Details).
#' @param p_criterion Information criterion used when `tweedie_p = "auto"`:
#'   `"BIC"` (default) or `"AIC"`. Both use the log-likelihood and effective
#'   degrees of freedom reported by [mgcv::logLik.gam()].
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
#' @param size_factor_method Method used to compute the size-factor offset.
#'   `"scran"` (default) uses pooling-based deconvolution size factors from
#'   [scran::computeSumFactors()]; `"library_size"` uses total counts scaled to
#'   a median of one; `"spanorm"` uses a SpaNorm fit, either as gene- and
#'   spot-specific size factors or as normalized counts (see `spanorm_output`
#'   and Details). For `"scran"` and `"library_size"`, an existing
#'   `colData(input)$sizeFactor` is used as-is. For `"spanorm"`, an existing
#'   `colData(input)$sizeFactor` is used by SpaNorm as its base library size.
#'   The SpaNorm fit is stored in the returned object's metadata with a
#'   fingerprint of the data, settings, and `seed`; later calls (for example
#'   [plot_feature()]) reuse it only when that fingerprint matches, and refit
#'   otherwise. SpaNorm subsamples locations at random, so set `seed` for
#'   reproducible results.
#' @param spanorm_output How the SpaNorm fit enters TweedSpot when
#'   `size_factor_method = "spanorm"`: `"size_factor"` (default) models the
#'   raw counts with the SpaNorm library-size (and batch) component as a gene-
#'   and spot-specific offset; `"counts"` models SpaNorm percentile-adjusted
#'   counts (PAC) with no further offset. Ignored for other methods.
#' @param spanorm_args Optional named list of additional arguments passed to
#'   `SpaNorm::SpaNorm()` when `size_factor_method = "spanorm"`, for example
#'   `list(sample.p = 0.25, df.tps = 6)`. `adj.method` is fixed to `"logpac"`.
#' @param padj_method Multiple-testing correction passed to [stats::p.adjust()].
#'   Default `"BY"`.
#' @param group Optional name of a discrete `colData(input)` column, such as a
#'   cell-type label. When given, each feature is tested for spatial variation
#'   *within* each group (see Details), and results are stored per feature and
#'   group instead of in the usual `rowData()` columns. Default `NULL`.
#' @param min_group_size Groups with fewer locations than this are not tested.
#'   Default `50`.
#' @param group_adjust Multiple-testing adjustment across features and groups:
#'   `"stagewise"` (default), `"global"`, or `"per_group"`. See Details.
#' @param alpha Target significance level for group results: used by the
#'   stage-wise adjustment and to summarise `rowData(input)$tweedspot_group_sig`.
#'   Default `0.05`.
#' @param seed Optional whole number. When given, the random-number generator
#'   is seeded with it for the duration of the call (and restored afterwards),
#'   so steps that use random numbers, such as SpaNorm's subsampling of
#'   locations, give identical results on every run. The seed is stored and
#'   reused by [plot_feature()] and [plot_top_features()]. Default `NULL`.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object controlling
#'   parallelization across genes. Default [BiocParallel::bpparam()], which
#'   uses the currently registered BiocParallel backend.
#' @param verbose Logical; if `TRUE`, print lightweight progress messages before
#'   and after gene-wise model fitting. Default `TRUE`.
#'
#' @details
#' **Normalization.** TweedSpot models raw counts with a log-link GAM and a
#' size-factor offset, so only the relative size factors across spots matter.
#' `scuttle::logNormCounts()` is not offered as a separate option because,
#' without user-supplied size factors, it uses library-size factors
#' (`colSums / mean(colSums)`) that differ from `"library_size"` only by a
#' global constant, which the model intercept absorbs; with existing size
#' factors it reuses them. Its log-transformed output is also incompatible with
#' modelling raw counts. `"spanorm"` (Salim et al., 2025, Genome Biology)
#' estimates a spatially smooth, gene-specific library-size effect jointly with
#' spatial biology, without domain or cluster labels. TweedSpot treats that
#' estimated effect as known, so any error in separating library size from
#' biology remains in the data as spatial structure. When library size is
#' itself spatially structured, this can make the spatial test strongly
#' anti-conservative (in simulations, most null genes were called
#' significant). Use `"spanorm"` as an exploratory option and check its
#' calibration on your data; `"scran"` remains the default. SpaNorm is specific
#' to spatial transcriptomics count data; `"scran"` and `"library_size"` are
#' the modality-agnostic options.
#'
#' SpaNorm fits, per gene, a negative binomial model whose log mean is
#' \eqn{\gamma_g + f_g(x, y) + \{a_g + h_g(x, y)\}\log LS + \text{batch}},
#' where \eqn{f_g} is spatial biology and the second term is a spatially
#' varying library-size effect. With `spanorm_output = "size_factor"`,
#' TweedSpot models raw counts and uses the library-size and batch part of
#' this predictor as the offset, \eqn{\log s_{gi}}, so the spatial biology is
#' tested by the TweedSpot smooth rather than removed or reused. With
#' `spanorm_output = "counts"`, TweedSpot models SpaNorm's percentile-adjusted
#' counts, \eqn{2^{\text{logpac}} - 1}, which are non-negative integers with
#' the library-size effect removed, using a unit offset. The size-factor form
#' is the default because it keeps the original count likelihood.
#'
#' **Automatic variance power.** With `tweedie_p = "auto"`, each gene is
#' fitted separately under several candidate variance functions
#' \eqn{V(\mu) = \phi\mu^p}: Gaussian (\eqn{p = 0}), Tweedie with \eqn{p}
#' estimated in \eqn{(1, 2)}, and, when the gene has no zeros, Gamma
#' (\eqn{p = 2}) and inverse Gaussian (\eqn{p = 3}). Every candidate shares the
#' log link, spatial smooth, covariates, and offset, and estimates its own
#' dispersion, so they differ only in the variance function. Poisson
#' (\eqn{p = 1}) is not a candidate: its dispersion is fixed at 1, and choosing
#' it on the same fit that is tested inflates false positives for genes whose
#' noise happens to look spatial. It remains available as `tweedie_p = 1`. Each candidate is scored by
#' \deqn{\mathrm{BIC} = -2\ell + \log(n)\,\mathrm{df}}
#' (or \eqn{\mathrm{AIC} = -2\ell + 2\,\mathrm{df}} with
#' `p_criterion = "AIC"`), where \eqn{\ell} is the maximized
#' log-likelihood and df the effective degrees of freedom from
#' [mgcv::logLik.gam()]: smooth EDF, plus one for an estimated scale, plus one
#' for an estimated Tweedie power. The candidate with the smallest value is used
#' for the spatial test, and its power is reported in
#' `rowData(input)$tweedspot_p` and `rowData(input)$tweedspot_family`. For
#' integer data the estimated Tweedie power is restricted to \eqn{p \ge 1.1},
#' because as \eqn{p \to 1} the Tweedie density concentrates on spikes at the
#' integers, which inflates its likelihood and the spatial test.
#'
#' **Group-specific testing.** With `group = "cell_type"`, the locations of each
#' group are fitted separately with the same model (`s(x, y)`, covariates,
#' offset, and family settings) to ask whether a feature varies in space among
#' locations of that group. A feature that differs between groups but is flat
#' within each group, a pure composition effect, is therefore not called
#' spatially variable. Size factors and the coordinate scaling are computed
#' once on all locations so that every group shares one scale; the spatial
#' basis size is capped by each group's number of locations; covariates that
#' are constant within a group are dropped for that group; and a feature is
#' tested in a group only if it is detected in at least as many of its
#' locations as the smooth has basis functions. Groups with fewer than
#' `min_group_size` locations are skipped.
#'
#' Group p-values are adjusted with `group_adjust`. `"stagewise"` first screens
#' each feature with the Cauchy combination of its group p-values and adjusts
#' those across features with `padj_method`, then Holm-adjusts the group
#' p-values within the feature at the reduced level `alpha * R / m` (R of m
#' features pass the screen); a feature-group pair's adjusted p-value is the
#' larger of the screening and confirmation values, so it is significant at
#' `alpha` only if the feature passes the screen and the group passes
#' confirmation. Stage-wise adjusted p-values are calibrated for the `alpha`
#' used in the run. `"global"` adjusts all pairs
#' together and `"per_group"` adjusts within each group. The full table, one
#' row per feature and group, is stored in
#' `S4Vectors::metadata(input)$tweedspot_group$results` and returned by
#' [get_tweedspot_results()] with its `group` argument; `rowData(input)` gains
#' `tweedspot_group_padj` (feature-level adjusted p-value) and
#' `tweedspot_group_sig` (groups significant at `alpha`).
#'
#' @return The input `SpatialExperiment` with per-gene statistic, p-value,
#'   adjusted p-value, smooth effective degrees of freedom, deviance
#'   explained, fitted variance power (`tweedspot_p`), and the corresponding
#'   family (`tweedspot_family`) written to [SummarizedExperiment::rowData()].
#'   The model settings are stored in `S4Vectors::metadata(input)$tweedspot`
#'   and reused by [plot_feature()] and [plot_top_features()].
#'
#' @references
#' Lun ATL, Bach K, Marioni JC (2016). Pooling across cells to normalize
#' single-cell RNA sequencing data with many zero counts. *Genome Biology*
#' 17:75.
#'
#' Salim A, Bhuva DD, Chen C, et al. (2025). SpaNorm: spatially-aware
#' normalization for spatial transcriptomics data. *Genome Biology*.
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
#' spe <- tweedspot(
#'   spe,
#'   family = "poisson",
#'   smooth_k = 4,
#'   size_factor_method = "library_size",
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#' SummarizedExperiment::rowData(spe)[, c("tweedspot_stat", "tweedspot_pval")]
#' @export
tweedspot <- function(input,
                      assay_name = "counts",
                      covariates = NULL,
                      two_part = FALSE,
                      family = "tw",
                      tweedie_p = NULL,
                      fit_method = "REML",
                      use_bam = FALSE,
                      bam_discrete = TRUE,
                      bam_nthreads = 1L,
                      smooth_k = NULL,
                      size_factor_method = c("scran", "library_size", "spanorm"),
                      spanorm_output = c("size_factor", "counts"),
                      spanorm_args = list(),
                      p_criterion = c("BIC", "AIC"),
                      padj_method = "BY",
                      group = NULL,
                      min_group_size = 50,
                      group_adjust = c("stagewise", "global", "per_group"),
                      alpha = 0.05,
                      seed = NULL,
                      BPPARAM = BiocParallel::bpparam(),
                      verbose = TRUE) {

  size_factor_method <- match.arg(size_factor_method)
  spanorm_output <- match.arg(spanorm_output)
  p_criterion <- match.arg(p_criterion)
  group_adjust <- match.arg(group_adjust)
  stopifnot(methods::is(input, "SpatialExperiment"))
  tweedspot_check_flag(verbose, "verbose")
  tweedspot_check_scalar(min_group_size, "min_group_size", min_value = 1)
  tweedspot_check_scalar(alpha, "alpha", min_value = 0, max_value = 1)
  restore_rng <- tweedspot_set_seed(seed)
  on.exit(restore_rng(), add = TRUE)
  labels <- if (is.null(group)) NULL else tweedspot_group_labels(input, group)
  # Start from a clean object so results from an earlier run (ordinary or
  # group-specific) can never be mixed with this run's results.
  input <- tweedspot_clear_results(input, verbose = verbose)

  family_arg <- family
  fit_input <- tweedspot_prepare_input(
    input = input,
    assay_name = assay_name,
    covariates = covariates,
    size_factor_method = size_factor_method,
    spanorm_args = spanorm_args,
    spanorm_output = spanorm_output,
    verbose = verbose,
    seed = seed
  )
  family <- tweedspot_resolve_family(
    family = family,
    tweedie_p = tweedie_p,
    Y = fit_input$Y
  )
  tweedspot_message_start(
    verbose = verbose,
    Y = fit_input$Y,
    BPPARAM = BPPARAM
  )
  if (is.null(labels)) {
    res <- tweedspot_agnostic(
      fit_input$Y,
      fit_input$coords,
      fit_input$libsz,
      fit_input$covariates,
      two_part,
      family,
      fit_method,
      use_bam,
      bam_discrete,
      bam_nthreads,
      smooth_k,
      BPPARAM,
      p_criterion = p_criterion
    )
    input <- tweedspot_store_results(input, res, padj_method = padj_method)
    n_fitted <- sum(!is.na(res$pval))
  } else {
    if (verbose && anyNA(labels)) {
      message(sum(is.na(labels)), " location(s) with a missing `", group,
              "` label are not tested.")
    }
    res <- tweedspot_group_fit(
      fit_input = fit_input,
      labels = labels,
      min_group_size = min_group_size,
      two_part = two_part,
      family = family,
      fit_method = fit_method,
      use_bam = use_bam,
      bam_discrete = bam_discrete,
      bam_nthreads = bam_nthreads,
      smooth_k = smooth_k,
      BPPARAM = BPPARAM,
      p_criterion = p_criterion,
      verbose = verbose
    )
    if (!any(res$tested)) {
      stop("No feature could be tested in any group. Check `min_group_size`, ",
           "`smooth_k`, and filtering thresholds.")
    }
    res <- tweedspot_group_adjust(res, group_adjust, padj_method, alpha)
    input <- tweedspot_store_group_results(input, res, group, group_adjust, alpha)
    n_fitted <- sum(res$tested)
  }
  if (!is.null(fit_input$spanorm_fit)) {
    S4Vectors::metadata(input)$SpaNorm <- fit_input$spanorm_fit
    S4Vectors::metadata(input)$SpaNorm_fingerprint <- fit_input$spanorm_fingerprint
  }
  # Record the model settings so plot_feature() and plot_top_features() can
  # refit each feature exactly as it was tested.
  S4Vectors::metadata(input)$tweedspot <- list(
    assay_name = assay_name,
    covariates = covariates,
    family = family_arg,
    tweedie_p = tweedie_p,
    p_criterion = p_criterion,
    fit_method = fit_method,
    use_bam = use_bam,
    bam_discrete = bam_discrete,
    smooth_k = smooth_k,
    size_factor_method = size_factor_method,
    spanorm_output = spanorm_output,
    spanorm_args = spanorm_args,
    two_part = two_part,
    group = group,
    min_group_size = min_group_size,
    group_adjust = group_adjust,
    seed = seed,
    size_factors = if (is.matrix(fit_input$libsz)) NULL else fit_input$libsz,
    features = rownames(input),
    locations = colnames(input)
  )
  if (verbose) {
    message("Completed TweedSpot fits for ", n_fitted,
            if (is.null(labels)) " feature(s)." else " feature-group pair(s).")
  }
  input
}
