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
#' spe <- tweedspot(spe, family = "poisson", smooth_k = 4, BPPARAM = BiocParallel::SerialParam())
#' SummarizedExperiment::rowData(spe)[, c("tweedspot_stat", "tweedspot_pval")]
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
  tweedspot_check_flag(verbose, "verbose")

  fit_input <- tweedspot_prepare_input(
    input = input,
    assay_name = assay_name,
    covariates = covariates
  )
  tweedspot_message_start(
    verbose = verbose,
    Y = fit_input$Y,
    BPPARAM = BPPARAM
  )
  res <- tweedspot_agnostic(
    fit_input$Y,
    fit_input$coords,
    fit_input$libsz,
    fit_input$covariates,
    two_part,
    combine,
                            family, fit_method, use_bam, bam_discrete,
                            bam_nthreads, smooth_k, BPPARAM)
  input <- tweedspot_store_results(input, res, padj_method = padj_method)
  tweedspot_message_end(verbose = verbose, pval = res$pval)
  input
}
