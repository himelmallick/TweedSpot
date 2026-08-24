# TweedSpot

TweedSpot is an R package for detecting spatially variable genes (SVGs) in
spatial omics data with spatial generalized additive models using a Tweedie
count family.

## Installation

TweedSpot is intended to be installed from GitHub.

```r
install.packages("remotes")
remotes::install_github("himelmallick/TweedSpot", dependencies = TRUE)
```

Because TweedSpot depends on Bioconductor infrastructure, install missing
Bioconductor packages first if needed:

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
  smooth_k = 20
)
```

## Main Functions

TweedSpot currently exposes two user-facing functions:

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

## Output

TweedSpot writes these per-gene fields to `rowData(input)`:

- `tweedspot_stat`
- `tweedspot_pval`
- `tweedspot_padj`
- `tweedspot_edf`
- `tweedspot_dev_expl`

## Status

The package is under active development. A vignette and broader test coverage
can be added on top of this installable package skeleton.

## License

MIT
