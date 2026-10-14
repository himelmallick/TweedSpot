#' TweedSpot: Self-Adaptive Spatial Tweedie Models for Spatially Variable Gene Detection
#'
#' TweedSpot provides spatially variable gene detection for spatial omics data
#' using Tweedie generalized additive models with spatial smoothing.
#'
#' The package currently exposes `filter_features()` for preprocessing and
#' `tweedspot()` for spatially variable gene detection with Tweedie generalized
#' additive models.
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom scran computeSumFactors
#' @importFrom MatrixGenerics rowSums colSums rowMeans colMeans
#' @keywords internal
#' @name TweedSpot-package
"_PACKAGE"

utils::globalVariables(c("value", "x", "y"))
