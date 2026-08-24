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
#'   p-values. Default `FALSE`.
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
#' @param verbose Logical; if `TRUE`, print lightweight progress messages before
#'   and after gene-wise model fitting. Default `TRUE`.
#'
#' @return The input `SpatialExperiment` with per-gene statistic, p-value,
#'   adjusted p-value, smooth effective degrees of freedom, and deviance
#'   explained written to [SummarizedExperiment::rowData()].
#'
#' @examples
#' \dontrun{
#' library(STexampleData)
#' spe <- ST_mouseOB()
#' spe <- tweedspot(spe)
#' head(SummarizedExperiment::rowData(spe))
#' }
#' @export
tweedspot <- function(input,
                      assay_name = "counts",
                      covariates = NULL,
                      two_part = FALSE,
                      combine = c("CCT", "stouffer", "min"),
                      family = "tw",
                      fit_method = "REML",
                      use_bam = FALSE,
                      bam_discrete = TRUE,
                      bam_nthreads = 1L,
                      smooth_k = NULL,
                      padj_method = "BY",
                      BPPARAM = BiocParallel::bpparam(),
                      verbose = TRUE) {

  combine <- match.arg(combine)
  stopifnot(methods::is(input, "SpatialExperiment"))

  Y <- as.matrix(SummarizedExperiment::assay(input, assay_name))
  coords <- scale(SpatialExperiment::spatialCoords(input))
  genes  <- rownames(input)
  libsz  <- tweedspot_libsize(input, Y)
  covariates <- tweedspot_covariates(input, covariates)

  if (ncol(coords) < 2) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be `TRUE` or `FALSE`.")
  }
  if (verbose) {
    message(
      "Running TweedSpot on ", nrow(Y), " genes across ", ncol(Y), " spatial locations ",
      "with ", BiocParallel::bpnworkers(BPPARAM), " worker(s)."
    )
  }
  res <- tweedspot_agnostic(Y, coords, libsz, covariates, two_part, combine,
                            family, fit_method, use_bam, bam_discrete,
                            bam_nthreads, smooth_k, BPPARAM)

  SummarizedExperiment::rowData(input)$tweedspot_stat <- res$stat
  SummarizedExperiment::rowData(input)$tweedspot_pval <- res$pval
  SummarizedExperiment::rowData(input)$tweedspot_padj <-
    stats::p.adjust(res$pval, method = padj_method)
  SummarizedExperiment::rowData(input)$tweedspot_edf <- res$edf
  SummarizedExperiment::rowData(input)$tweedspot_dev_expl <- res$dev_expl
  if (verbose) {
    n_sig <- sum(!is.na(res$pval))
    message("Completed TweedSpot fits for ", n_sig, " gene(s).")
  }
  input
}
