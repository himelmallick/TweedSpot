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

spe <- filter_features(
  spe,
  min_total = 5,
  min_pct_locations = 5,
  min_locations = 3,
  min_mean = 0.1,
  smooth_k = 20,                # same basis size as tweedspot() below
  exclude_mito = TRUE,          # transcriptomics only
  max_excluded_fraction = 0.2,  # drop locations with >20% mitochondrial signal
  outlier_nmads = 3,            # drop low-quality locations
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

- `filter_features()`: remove or flag low-information features (genes,
  proteins, peaks, metabolites), features with too few detected locations to
  support the spatial smooth, named controls, and low-quality locations;
  `action = "flag"` annotates instead of dropping. Mitochondrial and ribosomal
  gene removal is available for transcriptomics
- `get_tweedspot_results()`: extract and sort TweedSpot results from `rowData()`
- `top_spatial_genes()`: return the top-ranked spatially variable genes
- `tweedspot()`: run spatially variable gene detection on a
  `SpatialExperiment`
- `plot_feature()`: observed signal, TweedSpot fit, and fitted spatial effect
  for one feature, refitted with the settings of the `tweedspot()` run
- `plot_top_features()`: spatial effects of the `n_top` most significant
  features on a shared scale

Both plotting functions save PDF/PNG/SVG/TIFF files when given `output_dir`:

```r
plot_feature(spe_res, "GENE1", output_dir = "figures", formats = c("pdf", "png"))
plot_top_features(spe_res, n_top = 9, output_dir = "figures")
```

Key arguments:

- `seed`: optional whole number that makes runs reproducible; it is stored
  and reused by the plotting functions, and the caller's random-number state
  is restored afterwards

- `input`: a `SpatialExperiment`
- `assay_name`: assay to model, typically `"counts"`
- `covariates`: optional one-sided formula or character vector of `colData(input)` variables
- `two_part`: optional two-part agnostic model
- `use_bam`: use `mgcv::bam()` for faster fitting on larger datasets
- `smooth_k`: optional smooth basis size to trade flexibility for speed
- `size_factor_method`: `"scran"` (default) or `"library_size"` (total
  counts) offsets, both modality-agnostic, or `"spanorm"` for spatially-aware
  SpaNorm normalization of spatial transcriptomics counts
- `spanorm_output`: with `"spanorm"`, use SpaNorm as gene- and spot-specific
  size factors on raw counts (`"size_factor"`, default) or model SpaNorm
  percentile-adjusted counts (`"counts"`)
- `tweedie_p`: optionally select an `mgcv`-supported Tweedie variance power
  mapping, including estimated Tweedie (`NULL`), Gaussian (`0`), Poisson (`1`),
  fixed Tweedie powers in `(1, 2]`, or inverse Gaussian (`3`) for strictly
  positive responses

## Cell-Type-Specific Testing

For single-cell-resolution data (Xenium, CosMx, MERFISH, segmented imaging),
`group` tests each feature for spatial variation *within* each cell type, so a
marker that only reflects where a cell type sits is not called:

```r
spe <- filter_features(spe, group = "cell_type", smooth_k = 20)
spe <- tweedspot(spe, group = "cell_type", smooth_k = 20)
get_tweedspot_results(spe, group = "T cell")
```

Each type is fitted with its own spatial smooth on shared size factors and
coordinates; types smaller than `min_group_size` are skipped; and p-values are
adjusted stage-wise by default (`group_adjust = "stagewise"`, or `"global"` /
`"per_group"`).

## SpaNorm Normalization

For spatial transcriptomics, library size can itself vary with tissue
structure. [SpaNorm](https://bioconductor.org/packages/SpaNorm) estimates a
spatially smooth library-size effect separately from spatial biology, without
domain labels, and TweedSpot can use that fit as an optional alternative to
scran. TweedSpot treats SpaNorm's estimate as known, so when library size is
spatially structured the spatial test can be strongly anti-conservative; use
it as an exploratory option and check calibration on your data. TweedSpot can
use a SpaNorm fit in two ways:

```r
# Raw counts + SpaNorm library-size component as the offset
spe_res <- tweedspot(spe, size_factor_method = "spanorm")

# SpaNorm percentile-adjusted counts, no further offset
spe_res <- tweedspot(spe, size_factor_method = "spanorm",
                     spanorm_output = "counts")
```

SpaNorm subsamples spots when fitting, so pass `seed` (for example
`tweedspot(spe, size_factor_method = "spanorm", seed = 1)`) for reproducible
results. The fit is stored in `metadata(spe_res)$SpaNorm` with a fingerprint of
the data and settings, and reused by later calls such as `plot_feature()` only
when that fingerprint matches. Install SpaNorm with
`BiocManager::install("SpaNorm")`.

## Tweedie Index Parameter

TweedSpot uses `mgcv` as its modeling engine, so the supported Tweedie index
settings follow the families that can be represented cleanly inside
`mgcv::gam()` and `mgcv::bam()`. The `tweedie_p` argument controls this behavior
directly.

By default, `tweedie_p = NULL`, which fits an adaptive Tweedie GAM using
`mgcv::tw()`. In this mode, the Tweedie index is estimated during model fitting.
This is the recommended default for count-like spatial transcriptomics data
because it lets each gene adapt to its own mean-variance relationship.

Set `tweedie_p = "auto"` to choose the variance power per gene: TweedSpot fits
the Gaussian, estimated Tweedie, Gamma and inverse Gaussian models that are
valid for that gene's data, and keeps the one with the lowest BIC (or AIC via
`p_criterion = "AIC"`). Poisson is not a candidate because its fixed dispersion
makes the selected test anti-conservative; use `tweedie_p = 1` to fit it
explicitly. The choice is stored in
`rowData()$tweedspot_p` and `rowData()$tweedspot_family`.

You can also provide fixed values:

| `tweedie_p` | Model family used by TweedSpot | Typical use |
| --- | --- | --- |
| `NULL` | `mgcv::tw()` | Estimate the Tweedie index during fitting |
| `"auto"` | best of the candidates below by BIC/AIC | Let each gene pick its variance power |
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
