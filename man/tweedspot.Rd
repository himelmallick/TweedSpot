\name{tweedspot}
\alias{tweedspot}
\title{Detect spatially variable genes with self-adaptive spatial Tweedie models}
\usage{
tweedspot(
  input,
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
  BPPARAM = BiocParallel::bpparam()
)
}
\arguments{
  \item{input}{A \code{SpatialExperiment::SpatialExperiment} with a count assay and spatial coordinates.}

  \item{assay_name}{Name of the count assay to model. Default \code{"counts"}.}

  \item{W}{Optional spots-by-K matrix of cell-type proportions. Supplying \code{W} triggers proportion-based cell-type-specific detection.}

  \item{celltype}{Optional factor of length \code{ncol(input)} giving a cell-type label per location. Supplying \code{celltype} triggers label-based cell-type-specific detection. Ignored if \code{W} is supplied.}

  \item{two_part}{Logical; if \code{TRUE}, fit the two-part agnostic model and combine count and detection p-values.}

  \item{combine}{P-value combination rule for the two-part model: \code{"CCT"}, \code{"stouffer"}, or \code{"min"}.}

  \item{family}{\code{mgcv} family for the count model. Default \code{"tw"}.}

  \item{fit_method}{Smoothing parameter estimation method passed to \code{mgcv}, such as \code{"REML"} or \code{"GCV.Cp"}.}

  \item{use_bam}{Logical; use \code{mgcv::bam()} instead of \code{mgcv::gam()} for scalability.}

  \item{bam_discrete}{Logical; when \code{use_bam = TRUE}, use \code{mgcv::bam(discrete = TRUE)} for faster fitting.}

  \item{bam_nthreads}{Integer; number of threads to pass to \code{mgcv::bam()}.}

  \item{smooth_k}{Optional basis dimension for the spatial smooth. Smaller values can reduce runtime.}

  \item{padj_method}{Multiple-testing correction method passed to \code{stats::p.adjust()}.}

  \item{BPPARAM}{A \code{BiocParallel::BiocParallelParam} object controlling parallelization across genes. Defaults to the currently registered BiocParallel backend.}
}
\value{
The input \code{SpatialExperiment} with TweedSpot results added.

In agnostic mode, per-gene statistic, p-value, and adjusted p-value are written
to \code{rowData(input)}. In cell-type-specific modes, genes-by-cell-type
p-value and adjusted p-value matrices are written to \code{rowData(input)}, and
a long-format results table is attached to \code{metadata(input)$tweedspot}.
}
\description{
\code{tweedspot()} is the main entry point for TweedSpot. It fits spatial
generalized additive models gene-by-gene using a Tweedie count family by
default and supports agnostic, proportion-based, and label-based spatially
variable gene detection.
}
\details{
The package currently exposes a single user-facing function. In agnostic mode,
\code{tweedspot()} fits a model of the form \code{s(x, y)} per gene and tests
the spatial smooth. In cell-type-specific modes, it fits either separate smooths
weighted by a cell-type proportion matrix or factor-by smooths using cell-type
labels.
}
\examples{
\dontrun{
library(STexampleData)
spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]
spe <- tweedspot(spe, assay_name = "counts")
head(SummarizedExperiment::rowData(spe))
}
}
