# TweedSpot

TweedSpot is an R package for detecting spatially variable genes (SVGs) in
spatial omics data with spatial generalized additive models using a Tweedie
count family.

TweedSpot is designed for the broader spatial transcriptomics community, where
users often want methods that model count data directly, remain flexible across
different spatial expression patterns, and integrate cleanly with established
analysis workflows. We are submitting TweedSpot to Bioconductor because we
think it can be broadly useful as a community-facing method for reproducible
SVG discovery, benchmarking, and downstream biological interpretation.

## Installation

Once TweedSpot is available through Bioconductor, the release version can be
installed with `BiocManager`:

```r
install.packages("BiocManager")
BiocManager::install("TweedSpot")
```

The development version is currently available from GitHub:

```r
install.packages("remotes")
remotes::install_github("himelmallick/TweedSpot", dependencies = TRUE)
```

If you need to install core Bioconductor dependencies separately, you can use:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "SpatialExperiment",
  "SummarizedExperiment",
  "BiocParallel",
  "STexampleData"
))
```

## Quick Start

```r
library(TweedSpot)
library(SpatialExperiment)
library(STexampleData)

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

spe_res <- tweedspot(
  input = spe,
  assay_name = "counts"
)

head(SummarizedExperiment::rowData(spe_res)[, c(
  "tweedspot_stat",
  "tweedspot_pval",
  "tweedspot_padj",
  "tweedspot_edf",
  "tweedspot_dev_expl"
)])
```

## Example With Filtering

```r
library(TweedSpot)
library(SpatialExperiment)
library(STexampleData)

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

spe <- filter_genes(
  spe,
  filter_genes_ncounts = 5,
  filter_genes_pcspots = 5,
  filter_genes_nspots = 3,
  filter_genes_mean = 0.1,
  exclude_mito = TRUE,
  drop_zero_spots = TRUE,
  report = TRUE
)

spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  two_part = FALSE,
  family = "tw",
  fit_method = "REML",
  use_bam = TRUE,
  bam_discrete = TRUE,
  size_factor_method = "scran",
  smooth_k = 20
)
```

## Main Functions

TweedSpot currently exposes the following user-facing functions:

- `filter_genes()`: prefilter low-information genes before model fitting
- `get_tweedspot_results()`: extract and sort TweedSpot results from `rowData()`
- `top_spatial_genes()`: return the top-ranked spatially variable genes
- `tweedspot()`: run spatially variable gene detection on a
  `SpatialExperiment`

Key arguments:

- `input`: a `SpatialExperiment`
- `assay_name`: assay to model, typically `"counts"`
- `covariates`: optional one-sided formula or character vector of `colData(input)` variables
- `two_part`: optional two-part agnostic model
- `use_bam`: use `mgcv::bam()` for faster fitting on larger datasets
- `smooth_k`: optional smooth basis size to trade flexibility for speed
- `size_factor_method`: use `scran` or total-count library-size offsets when
  `colData(input)$sizeFactor` is not already present
- `tweedie_p`: optionally select an `mgcv`-supported Tweedie variance power
  mapping, including estimated Tweedie (`NULL`), Gaussian (`0`), Poisson (`1`),
  fixed Tweedie powers in `(1, 2]`, or inverse Gaussian (`3`) for strictly
  positive responses

## Tweedie Index Parameter

TweedSpot uses `mgcv` as its modeling engine, so the supported Tweedie index
settings follow the families that can be represented cleanly inside
`mgcv::gam()` and `mgcv::bam()`. The `tweedie_p` argument controls this behavior
directly.

By default, `tweedie_p = NULL`, which fits an adaptive Tweedie GAM using
`mgcv::tw()`. In this mode, the Tweedie index is estimated during model fitting.
This is the recommended default for count-like spatial transcriptomics data
because it lets each gene adapt to its own mean-variance relationship.

You can also provide fixed values:

| `tweedie_p` | Model family used by TweedSpot | Typical use |
| --- | --- | --- |
| `NULL` | `mgcv::tw()` | Estimate the Tweedie index during fitting |
| `0` | `stats::gaussian()` | Continuous approximately normal outcomes |
| `1` | `stats::poisson()` | Poisson count model |
| `1 < p < 2` | `mgcv::Tweedie(p = p)` | Compound Poisson-Gamma Tweedie model with a point mass at zero |
| `2` | `mgcv::Tweedie(p = 2)` | Gamma boundary case |
| `3` | `stats::inverse.gaussian()` | Strictly positive continuous outcomes |

Values between `0` and `1`, values between `2` and `3`, negative values, and
values greater than `3` are not supported in the current `mgcv` implementation
used by TweedSpot. The inverse Gaussian option (`tweedie_p = 3`) requires
strictly positive responses, so it is usually not appropriate for raw spatial
count matrices containing zeros.

Examples:

```r
# Estimate the Tweedie index per gene
spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  tweedie_p = NULL
)

# Fixed compound Poisson-Gamma Tweedie index
spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  tweedie_p = 1.5
)

# Poisson spatial GAM
spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  tweedie_p = 1
)

# Gaussian spatial GAM for continuous assays
spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  tweedie_p = 0
)
```

The older `family` argument is still available for users who want to pass an
`mgcv` family directly. When `tweedie_p` is supplied, it takes precedence over
`family`.

## Output

TweedSpot writes these per-gene fields to `rowData(input)`:

- `tweedspot_stat`
- `tweedspot_pval`
- `tweedspot_padj`
- `tweedspot_edf`
- `tweedspot_dev_expl`

## Package Notes

The package is built around the `SpatialExperiment` container and writes
gene-level inference back into `rowData()`. See the vignette for a fuller
workflow including filtering, parallel model fitting, result extraction, and
visualization.

## License

MIT
